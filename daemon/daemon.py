#!/usr/bin/env python3
"""
Dark Gene Pipeline Portal Daemon

Portal과 연동하여:
1. 분석 완료된 샘플 자동 업로드
2. 실시간 분석 상태 모니터링
3. Portal API를 통한 상태 보고
"""

import os
import sys
import time
import json
import logging
import requests
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import hashlib
import schedule
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

# Configuration
BASE_DIR = os.environ.get('PIPELINE_BASE_DIR', '/pipeline')
FASTQ_BASE = os.path.join(BASE_DIR, 'fastq')
ANALYSIS_BASE = os.path.join(BASE_DIR, 'analysis')
OUTPUT_BASE = os.path.join(BASE_DIR, 'output')
LOG_BASE = os.path.join(BASE_DIR, 'log')

PORTAL_URL = os.environ.get('PORTAL_URL', 'https://portal.genolyx.com')
API_KEY = os.environ.get('PORTAL_API_KEY', '')
INSTITUTION_ID = os.environ.get('INSTITUTION_ID', 'default')

# Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/var/log/daemon.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('DarkGeneDaemon')


class PortalAPI:
    """Portal API 클라이언트"""
    
    def __init__(self, base_url: str, api_key: str):
        self.base_url = base_url.rstrip('/')
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'Authorization': f'Bearer {api_key}',
            'Content-Type': 'application/json'
        })
    
    def get_analysis_summary(self) -> Dict:
        """Portal에서 분석 상태 조회"""
        try:
            response = self.session.get(f'{self.base_url}/api/v1/analysis/summary')
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Failed to get analysis summary: {e}")
            return {}
    
    def create_analysis_order(self, work_dir: str, sample_name: str) -> Optional[str]:
        """새 분석 주문 생성"""
        try:
            data = {
                'work_dir': work_dir,
                'sample_name': sample_name,
                'status': 'WAITING',
                'institution_id': INSTITUTION_ID,
                'created_at': datetime.now().isoformat()
            }
            response = self.session.post(f'{self.base_url}/api/v1/analysis/order', json=data)
            response.raise_for_status()
            result = response.json()
            return result.get('order_id')
        except Exception as e:
            logger.error(f"Failed to create order: {e}")
            return None
    
    def update_analysis_status(self, order_id: str, status: str, progress: int = 0, details: str = ''):
        """분석 상태 업데이트"""
        try:
            data = {
                'status': status,
                'progress': progress,
                'details': details,
                'updated_at': datetime.now().isoformat()
            }
            response = self.session.put(f'{self.base_url}/api/v1/analysis/order/{order_id}', json=data)
            response.raise_for_status()
            logger.info(f"Updated order {order_id}: {status} ({progress}%)")
            return True
        except Exception as e:
            logger.error(f"Failed to update status: {e}")
            return False
    
    def upload_file(self, order_id: str, file_path: str, file_type: str) -> bool:
        """파일 업로드"""
        try:
            if not os.path.exists(file_path):
                logger.warning(f"File not found: {file_path}")
                return False
            
            # 파일 메타데이터
            file_name = os.path.basename(file_path)
            file_size = os.path.getsize(file_path)
            
            # Multipart upload
            with open(file_path, 'rb') as f:
                files = {'file': (file_name, f, 'application/octet-stream')}
                data = {
                    'order_id': order_id,
                    'file_type': file_type,
                    'file_size': file_size
                }
                
                # Temporarily remove Content-Type for multipart
                headers = {'Authorization': f'Bearer {self.api_key}'}
                response = requests.post(
                    f'{self.base_url}/api/v1/analysis/upload',
                    files=files,
                    data=data,
                    headers=headers,
                    timeout=300
                )
                response.raise_for_status()
            
            logger.info(f"Uploaded {file_name} for order {order_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to upload file {file_path}: {e}")
            return False
    
    def complete_analysis(self, order_id: str, duration: str):
        """분석 완료 보고"""
        try:
            data = {
                'status': 'COMPLETED',
                'duration': duration,
                'completed_at': datetime.now().isoformat()
            }
            response = self.session.post(f'{self.base_url}/api/v1/analysis/order/{order_id}/complete', json=data)
            response.raise_for_status()
            logger.info(f"Completed order {order_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to complete order: {e}")
            return False


class AnalysisMonitor:
    """분석 상태 모니터링"""
    
    def __init__(self, portal_api: PortalAPI):
        self.portal_api = portal_api
        self.state_file = '/var/lib/daemon/state.json'
        self.state = self.load_state()
    
    def load_state(self) -> Dict:
        """상태 파일 로드"""
        if os.path.exists(self.state_file):
            try:
                with open(self.state_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"WARNING: state file read error: {e}", file=sys.stderr)
        return {'orders': {}}
    
    def save_state(self):
        """상태 파일 저장"""
        os.makedirs(os.path.dirname(self.state_file), exist_ok=True)
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2)
    
    def scan_samples(self) -> List[Dict]:
        """샘플 스캔"""
        samples = []
        
        if not os.path.exists(FASTQ_BASE):
            return samples
        
        for work_dir in os.listdir(FASTQ_BASE):
            work_path = os.path.join(FASTQ_BASE, work_dir)
            if not os.path.isdir(work_path):
                continue
            
            for sample_name in os.listdir(work_path):
                sample_path = os.path.join(work_path, sample_name)
                if not os.path.isdir(sample_path):
                    continue
                
                # 상태 확인
                completed_marker = os.path.join(sample_path, 'analysis.completed')
                is_completed = os.path.exists(completed_marker)
                
                # R1/R2 확인
                files = os.listdir(sample_path)
                has_fastq = any(f.endswith(('.fq.gz', '.fastq.gz')) for f in files)
                
                # Analysis 디렉토리 확인
                analysis_path = os.path.join(ANALYSIS_BASE, work_dir, sample_name)
                is_running = os.path.exists(analysis_path) and not is_completed
                
                # Output 디렉토리 확인
                output_path = os.path.join(OUTPUT_BASE, work_dir, sample_name)
                has_output = os.path.exists(output_path)
                
                samples.append({
                    'work_dir': work_dir,
                    'sample_name': sample_name,
                    'sample_path': sample_path,
                    'analysis_path': analysis_path,
                    'output_path': output_path,
                    'is_completed': is_completed,
                    'is_running': is_running,
                    'has_fastq': has_fastq,
                    'has_output': has_output
                })
        
        return samples
    
    def get_analysis_progress(self, work_dir: str, sample_name: str) -> int:
        """분석 진행률 계산"""
        trace_file = os.path.join(ANALYSIS_BASE, work_dir, sample_name, 'pipeline_info', 'trace.txt')
        
        if not os.path.exists(trace_file):
            return 0
        
        STAGES = {
            'ALIGN_AND_SORT': 25,
            'MARK_DUPLICATES': 30,
            'CALL_VARIANTS': 50,
            'EXPANSION_HUNTER': 70,
            'GCNV_COHORT_RUN': 80,
            'POSTPROCESS_GCNV': 85,
            'GENERATE_VISUAL_EVIDENCE': 90,
            'GENERATE_SUMMARY_REPORT': 95
        }
        
        progress = 0
        try:
            with open(trace_file, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    
                    full_name = parts[1]
                    process_name = full_name.split(' (')[0].split(':')[-1]
                    
                    pct = STAGES.get(process_name, 0)
                    if pct > progress:
                        progress = pct
        except (IOError, OSError, ValueError) as e:
            print(f"WARNING: trace parse error: {e}", file=sys.stderr)
        
        return progress
    
    def upload_results(self, order_id: str, work_dir: str, sample_name: str) -> bool:
        """결과 파일 업로드"""
        output_path = os.path.join(OUTPUT_BASE, work_dir, sample_name)
        
        if not os.path.exists(output_path):
            logger.warning(f"Output path not found: {output_path}")
            return False
        
        upload_map = {
            'summary': ['*_summary_report.txt', '*_detailed_report.txt'],
            'snapshots': ['*_visual_report.html'],
            # *_annotated.vcf.gz : VEP CSQ 포함 (gene/HGVSc/effect 표시용) — 현재 파이프라인 출력 파일명
            # *_filtered.vcf.gz  : annotation 전 필터링 VCF (CSQ 없음, fallback용으로 유지)
            'vcf': [
                '*_annotated.vcf.gz', '*_annotated.vcf.gz.tbi',
                '*_filtered.vcf.gz', '*_filtered.vcf.gz.tbi',
            ],
            'bam': ['*.md.bam', '*.md.bam.bai'],
            'cnv': ['*_cnv.vcf.gz'],
            'sv': ['*_manta.vcf.gz'],
            'repeat': ['*_eh.vcf', '*_eh.json', '*.svg']
        }
        
        success_count = 0
        total_count = 0
        
        for file_type, patterns in upload_map.items():
            type_dir = os.path.join(output_path, file_type)
            if not os.path.exists(type_dir):
                continue
            
            for pattern in patterns:
                import glob
                for file_path in glob.glob(os.path.join(type_dir, pattern)):
                    total_count += 1
                    if self.portal_api.upload_file(order_id, file_path, file_type):
                        success_count += 1
        
        logger.info(f"Uploaded {success_count}/{total_count} files for order {order_id}")
        return success_count > 0
    
    def process_completed_samples(self):
        """완료된 샘플 처리"""
        samples = self.scan_samples()
        
        for sample in samples:
            if not sample['is_completed'] or not sample['has_output']:
                continue
            
            sample_key = f"{sample['work_dir']}/{sample['sample_name']}"
            
            # 이미 업로드한 샘플은 스킵
            if sample_key in self.state['orders']:
                order_info = self.state['orders'][sample_key]
                if order_info.get('uploaded', False):
                    continue
            
            # Portal에 order 생성 (또는 조회)
            if sample_key not in self.state['orders']:
                order_id = self.portal_api.create_analysis_order(
                    sample['work_dir'],
                    sample['sample_name']
                )
                if not order_id:
                    continue
                
                self.state['orders'][sample_key] = {
                    'order_id': order_id,
                    'created_at': datetime.now().isoformat(),
                    'uploaded': False
                }
                self.save_state()
            
            order_info = self.state['orders'][sample_key]
            order_id = order_info['order_id']
            
            # 상태 업데이트: 완료
            self.portal_api.update_analysis_status(order_id, 'COMPLETED', 100)
            
            # 결과 파일 업로드
            logger.info(f"Uploading results for {sample_key}")
            if self.upload_results(order_id, sample['work_dir'], sample['sample_name']):
                # Duration 계산
                log_file = os.path.join(LOG_BASE, sample['work_dir'], sample['sample_name'], 'nextflow.log')
                duration = self.calculate_duration(log_file)
                
                # 완료 보고
                self.portal_api.complete_analysis(order_id, duration)
                
                # 상태 저장
                order_info['uploaded'] = True
                order_info['uploaded_at'] = datetime.now().isoformat()
                self.save_state()
    
    def process_running_samples(self):
        """실행 중인 샘플 처리"""
        samples = self.scan_samples()
        
        for sample in samples:
            if not sample['is_running']:
                continue
            
            sample_key = f"{sample['work_dir']}/{sample['sample_name']}"
            
            # Order 생성 (없으면)
            if sample_key not in self.state['orders']:
                order_id = self.portal_api.create_analysis_order(
                    sample['work_dir'],
                    sample['sample_name']
                )
                if not order_id:
                    continue
                
                self.state['orders'][sample_key] = {
                    'order_id': order_id,
                    'created_at': datetime.now().isoformat(),
                    'uploaded': False
                }
                self.save_state()
            
            order_info = self.state['orders'][sample_key]
            order_id = order_info['order_id']
            
            # 진행률 업데이트
            progress = self.get_analysis_progress(sample['work_dir'], sample['sample_name'])
            self.portal_api.update_analysis_status(order_id, 'RUNNING', progress)
    
    def calculate_duration(self, log_file: str) -> str:
        """분석 시간 계산"""
        if not os.path.exists(log_file):
            return "Unknown"
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
                # Nextflow duration 파싱
                import re
                match = re.search(r'Duration\s*:\s*(.+)', content)
                if match:
                    return match.group(1).strip()
        except (IOError, OSError) as e:
            print(f"WARNING: log parse error: {e}", file=sys.stderr)
        
        return "Unknown"
    
    def get_summary(self) -> Dict:
        """분석 상태 요약"""
        samples = self.scan_samples()
        
        today = datetime.now().date().isoformat()
        
        summary = {
            'today_requested': 0,
            'running': 0,
            'queue_waiting': 0,
            'today_completed': 0,
            'today_failed': 0,
            'total_requested': len(self.state['orders']),
            'total_completed': 0,
            'total_failed': 0
        }
        
        # 실시간 스캔
        for sample in samples:
            sample_key = f"{sample['work_dir']}/{sample['sample_name']}"
            
            if sample['is_running']:
                summary['running'] += 1
            elif sample['has_fastq'] and not sample['is_completed']:
                summary['queue_waiting'] += 1
            elif sample['is_completed']:
                summary['today_completed'] += 1
        
        # State 기반 카운트
        for sample_key, order_info in self.state['orders'].items():
            if order_info.get('uploaded', False):
                summary['total_completed'] += 1
                
                # 오늘 완료 여부
                uploaded_at = order_info.get('uploaded_at', '')
                if uploaded_at.startswith(today):
                    summary['today_completed'] += 1
        
        return summary


class CompletionWatcher(FileSystemEventHandler):
    """분석 완료 감지 (실시간)"""
    
    def __init__(self, monitor: AnalysisMonitor):
        self.monitor = monitor
    
    def on_created(self, event):
        """파일 생성 이벤트"""
        if event.is_directory:
            return
        
        if event.src_path.endswith('analysis.completed'):
            logger.info(f"Detected completion: {event.src_path}")
            # 즉시 처리
            self.monitor.process_completed_samples()


def main():
    """메인 루프"""
    logger.info("=== Dark Gene Pipeline Daemon Starting ===")
    
    # Portal API 초기화
    if not API_KEY:
        logger.error("PORTAL_API_KEY not set!")
        sys.exit(1)
    
    portal_api = PortalAPI(PORTAL_URL, API_KEY)
    monitor = AnalysisMonitor(portal_api)
    
    # File system watcher 설정
    observer = Observer()
    handler = CompletionWatcher(monitor)
    
    if os.path.exists(FASTQ_BASE):
        observer.schedule(handler, FASTQ_BASE, recursive=True)
        observer.start()
        logger.info(f"Watching {FASTQ_BASE} for completion markers")
    
    # 주기적 작업 스케줄
    schedule.every(30).seconds.do(monitor.process_completed_samples)
    schedule.every(10).seconds.do(monitor.process_running_samples)
    
    logger.info("Daemon started successfully")
    
    try:
        while True:
            schedule.run_pending()
            time.sleep(1)
    except KeyboardInterrupt:
        logger.info("Shutting down...")
        observer.stop()
        observer.join()
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
