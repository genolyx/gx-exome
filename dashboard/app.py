import os
import subprocess
import sys
import threading
import time
import shutil
import psutil
from datetime import datetime
from flask import Flask, render_template, request, jsonify, send_from_directory, Response, redirect

app = Flask(__name__)

# --- Configuration ---
BASE_DIR = os.getenv('BASE_DIR', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
FASTQ_BASE = os.path.join(BASE_DIR, 'fastq')
ANALYSIS_BASE = os.path.join(BASE_DIR, 'analysis')
OUTPUT_BASE = os.path.join(BASE_DIR, 'output')
LOG_BASE = os.path.join(BASE_DIR, 'log')
DAEMON_API_URL = os.getenv('DAEMON_API_URL', 'http://localhost:8080')

app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024 * 1024

# Ensure directories exist
os.makedirs(FASTQ_BASE, exist_ok=True)
os.makedirs(ANALYSIS_BASE, exist_ok=True)
os.makedirs(OUTPUT_BASE, exist_ok=True)
os.makedirs(LOG_BASE, exist_ok=True)

# Global State
pipeline_process = None
pipeline_lock = threading.Lock()
current_analysis = {}  # {work_dir, sample_name, log_file}

def get_current_workdir():
    """Get current work directory in YYMM format."""
    now = datetime.now()
    return now.strftime("%y%m")

def list_fastq_samples(show_completed=False):
    """
    List all sample directories in fastq/<work_dir>/<sample_name>/
    Returns: [(work_dir, sample_name, is_completed, has_r1r2), ...]
    """
    samples = []
    
    if not os.path.exists(FASTQ_BASE):
        return samples
    
    # Iterate work_dirs (e.g., 2601, 2602)
    for work_dir in sorted(os.listdir(FASTQ_BASE)):
        work_path = os.path.join(FASTQ_BASE, work_dir)
        if not os.path.isdir(work_path):
            continue
        
        # Iterate sample_name dirs
        for sample_name in sorted(os.listdir(work_path)):
            sample_path = os.path.join(work_path, sample_name)
            if not os.path.isdir(sample_path):
                continue
            
            # Check completion marker
            completed_marker = os.path.join(sample_path, 'analysis.completed')
            is_completed = os.path.exists(completed_marker)
            
            # Skip completed if show_completed=False
            if is_completed and not show_completed:
                continue
            
            # Check R1/R2 existence
            files = os.listdir(sample_path)
            r1_exists = any('_R1' in f or '_1.fq' in f or '_1.fastq' in f for f in files)
            r2_exists = any('_R2' in f or '_2.fq' in f or '_2.fastq' in f for f in files)
            has_r1r2 = r1_exists and r2_exists
            
            samples.append({
                'work_dir': work_dir,
                'sample_name': sample_name,
                'is_completed': is_completed,
                'has_r1r2': has_r1r2,
                'path': f'{work_dir}/{sample_name}'
            })
    
    return samples

def get_fastq_pairs(work_dir, sample_name):
    """
    Get R1/R2 file pairs from fastq/<work_dir>/<sample_name>/
    Returns: (r1_path, r2_path) or (None, None) if not found
    """
    sample_path = os.path.join(FASTQ_BASE, work_dir, sample_name)
    if not os.path.exists(sample_path):
        return None, None
    
    files = [f for f in os.listdir(sample_path) if f.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fastq'))]
    
    r1_files = [f for f in files if '_R1' in f or '_1.fq' in f or '_1.fastq' in f]
    r2_files = [f for f in files if '_R2' in f or '_2.fq' in f or '_2.fastq' in f]
    
    if len(r1_files) == 1 and len(r2_files) == 1:
        return os.path.join(sample_path, r1_files[0]), os.path.join(sample_path, r2_files[0])
    
    return None, None

def get_pipeline_status():
    """Get current pipeline status."""
    global pipeline_process
    if pipeline_process is None:
        return "IDLE"
    
    poll = pipeline_process.poll()
    if poll is None:
        return "RUNNING"
    elif poll == 0:
        return "SUCCESS"
    else:
        return "FAILED"

def check_active_pipeline():
    """Check for active Nextflow/BWA processes."""
    global pipeline_process
    
    if pipeline_process is not None and pipeline_process.poll() is None:
        return True, "RUNNING"

    try:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                cmd = proc.info['cmdline']
                if cmd:
                    cmd_str = ' '.join(cmd)
                    if 'nextflow' in cmd_str or 'bwa' in cmd_str or 'gatk' in cmd_str:
                        if 'app.py' not in cmd_str and 'grep' not in cmd_str:
                            return True, "RUNNING (Attached)"
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
    except Exception:
        pass
        
    return False, "IDLE"

@app.route('/')
def index():
    """Main page."""
    # List samples (only incomplete by default)
    samples = list_fastq_samples(show_completed=False)
    
    # List results from output/
    results = []
    if os.path.exists(OUTPUT_BASE):
        for work_dir in os.listdir(OUTPUT_BASE):
            work_path = os.path.join(OUTPUT_BASE, work_dir)
            if not os.path.isdir(work_path):
                continue
            
            for sample_name in os.listdir(work_path):
                sample_path = os.path.join(work_path, sample_name)
                if not os.path.isdir(sample_path):
                    continue
                
                # Check for summary/snapshots
                summary_dir = os.path.join(sample_path, 'summary')
                snapshot_dir = os.path.join(sample_path, 'snapshots')
                
                if os.path.exists(summary_dir):
                    for f in os.listdir(summary_dir):
                        if f.endswith(('.txt', '.html')):
                            results.append({
                                'name': f,
                                'type': 'Report',
                                'path': f'{work_dir}/{sample_name}/summary',
                                'sample': f'{work_dir}/{sample_name}'
                            })
                
                if os.path.exists(snapshot_dir):
                    for f in os.listdir(snapshot_dir):
                        if f.endswith(('.html', '.png', '.svg')):
                            results.append({
                                'name': f,
                                'type': 'Snapshot',
                                'path': f'{work_dir}/{sample_name}/snapshots',
                                'sample': f'{work_dir}/{sample_name}'
                            })
    
    return render_template('index.html', samples=samples, results=results, status=get_pipeline_status())

@app.route('/list_samples', methods=['POST'])
def list_samples():
    """List samples (with show_completed toggle)."""
    data = request.json
    show_completed = data.get('show_completed', False)
    
    samples = list_fastq_samples(show_completed=show_completed)
    return jsonify({'samples': samples})

@app.route('/create_workdir', methods=['POST'])
def create_workdir():
    """Create a new work directory (YYMM format)."""
    data = request.json
    work_dir = data.get('work_dir', get_current_workdir())
    
    # Create all base directories
    for base in [FASTQ_BASE, ANALYSIS_BASE, OUTPUT_BASE, LOG_BASE]:
        path = os.path.join(base, work_dir)
        os.makedirs(path, exist_ok=True)
    
    return jsonify({'message': f'Work directory {work_dir} created', 'work_dir': work_dir})

@app.route('/upload', methods=['POST'])
def upload_file():
    """Upload FASTQ files to specific sample directory."""
    if 'files[]' not in request.files:
        return jsonify({'error': 'No file part'}), 400
    
    work_dir = request.form.get('work_dir', get_current_workdir())
    sample_name = request.form.get('sample_name')
    
    if not sample_name:
        return jsonify({'error': 'No sample_name provided'}), 400
    
    # Create sample directory
    sample_path = os.path.join(FASTQ_BASE, work_dir, sample_name)
    os.makedirs(sample_path, exist_ok=True)
    
    files = request.files.getlist('files[]')
    saved_files = []
    
    for file in files:
        if file.filename == '':
            continue
        if file:
            filename = file.filename
            filepath = os.path.join(sample_path, filename)
            file.save(filepath)
            saved_files.append(filename)
    
    return jsonify({'message': f'Uploaded {len(saved_files)} files to {work_dir}/{sample_name}', 'files': saved_files})

@app.route('/start', methods=['POST'])
def start_pipeline():
    """Start analysis for selected sample directory."""
    global pipeline_process, current_analysis
    
    data = request.json
    work_dir = data.get('work_dir')
    sample_name = data.get('sample_name')
    
    if not work_dir or not sample_name:
        return jsonify({'error': 'work_dir and sample_name required!'}), 400

    with pipeline_lock:
        if get_pipeline_status() == "RUNNING":
            return jsonify({'error': 'Pipeline is already running!'}), 409
        
        # Verify R1/R2 existence
        r1, r2 = get_fastq_pairs(work_dir, sample_name)
        if not r1 or not r2:
            return jsonify({'error': f'R1/R2 pair not found in {work_dir}/{sample_name}'}), 400
        
        # Check if already completed
        fastq_path = os.path.join(FASTQ_BASE, work_dir, sample_name)
        completed_marker = os.path.join(fastq_path, 'analysis.completed')
        if os.path.exists(completed_marker):
            # Remove marker to allow re-analysis
            os.remove(completed_marker)
        
        # Setup directories
        analysis_dir = os.path.join(ANALYSIS_BASE, work_dir, sample_name)
        output_dir = os.path.join(OUTPUT_BASE, work_dir, sample_name)
        log_dir = os.path.join(LOG_BASE, work_dir, sample_name)
        
        os.makedirs(analysis_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        
        # Log file
        log_file = os.path.join(log_dir, 'nextflow.log')
        
        # Store current analysis info
        current_analysis = {
            'work_dir': work_dir,
            'sample_name': sample_name,
            'log_file': log_file,
            'fastq_path': fastq_path,
            'analysis_dir': analysis_dir,
            'output_dir': output_dir
        }
        
        # Run Pipeline (Background Thread)
        def run_and_finalize():
            global pipeline_process, current_analysis
            
            # Analysis log file (모든 로그를 여기에 저장)
            analysis_log = os.path.join(analysis_dir, f'{sample_name}.analysis.log')
            
            cmd = [
                'nextflow', 'run', '/app/bin/main.nf',
                '-ansi-log', 'false',
                '--fastq_dir', fastq_path,
                '--outdir', analysis_dir,
                '--output_dir', output_dir,
                '--sample_name', sample_name,
                '--gcnv_model', '/app/models/gcnv_12sample_cohort',
                '-work-dir', os.path.join(analysis_dir, 'work'),
                '-with-report', os.path.join(log_dir, 'report.html'),
                '-with-trace', os.path.join(log_dir, 'trace.txt'),
                '-with-timeline', os.path.join(log_dir, 'timeline.html')
            ]
            # Plain logs: one Submitted/Completed line per task; avoid NF 25 live progress table spam in files
            run_env = os.environ.copy()
            run_env['NXF_ANSI_LOG'] = 'false'
            run_env['NXF_ANSI_SUMMARY'] = 'false'
            run_env.setdefault('NO_COLOR', '1')

            try:
                # analysis.log에 모든 출력 저장
                with open(analysis_log, 'w') as f_analysis:
                    # nextflow.log에도 저장 (기존 호환성)
                    with open(log_file, 'w') as f_log:
                        pipeline_process = subprocess.Popen(
                        cmd,
                        cwd='/app/bin',
                        env=run_env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        bufsize=1
                        )
                                                # 실시간으로 두 파일에 모두 쓰기
                        for line in pipeline_process.stdout:
                            f_analysis.write(line)
                            f_analysis.flush()
                            f_log.write(line)
                            f_log.flush()
                        
                        pipeline_process.wait()
                
                # Post-Run Logic
                if pipeline_process.returncode == 0:
                    print(f"Pipeline Success for {work_dir}/{sample_name}")
                    # Mark as completed
                    with open(os.path.join(fastq_path, 'analysis.completed'), 'w') as f:
                        f.write(f"Completed at: {datetime.now().isoformat()}\n")
                else:
                    print(f"Pipeline Failed for {work_dir}/{sample_name}")
                        
            except Exception as e:
                print(f"Critical Pipeline Error: {e}")
                # 에러도 로그에 기록
                with open(analysis_log, 'a') as f:
                    f.write(f"\nERROR: {e}\n")
            finally:
                pipeline_process = None

        # Start thread
        thread = threading.Thread(target=run_and_finalize)
        thread.start()
        
        return jsonify({
            'message': f'Started analysis for {work_dir}/{sample_name}',
            'status': 'RUNNING',
            'work_dir': work_dir,
            'sample_name': sample_name
        })

@app.route('/stop', methods=['POST'])
def stop_pipeline():
    """Stop running pipeline."""
    global pipeline_process, current_analysis
    
    if pipeline_process and pipeline_process.poll() is None:
        try:
            parent = psutil.Process(pipeline_process.pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
            pipeline_process = None
        except Exception as e:
            print(f"Error killing process: {e}")

    # Kill external processes
    try:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                cmd_str = ' '.join(proc.info['cmdline'] or [])
                if 'nextflow' in cmd_str and 'main.nf' in cmd_str:
                    if proc.pid != os.getpid():
                        for child in proc.children(recursive=True):
                            child.kill()
                        proc.kill()
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                pass
    except Exception as e:
        print(f"WARNING: stop_pipeline error: {e}", file=sys.stderr)
    
    current_analysis = {}
    return jsonify({'message': 'Pipeline stopped.', 'status': 'IDLE'})

@app.route('/status')
def status():
    """Get pipeline status with progress."""
    global current_analysis
    
    is_running, state_str = check_active_pipeline()
    state = "RUNNING" if is_running else get_pipeline_status()
    
    # Check log for completion
    if state == "IDLE" and current_analysis.get('log_file'):
        log_file = current_analysis['log_file']
        if os.path.exists(log_file):
            try:
                with open(log_file, 'r') as f:
                    f.seek(0, os.SEEK_END)
                    f_size = f.tell()
                    f.seek(max(f_size - 1024, 0))
                    last_lines = f.read()
                    
                    if "Execution complete -- Goodbye" in last_lines:
                        state = "SUCCESS"
                    elif "ERROR" in last_lines:
                         state = "FAILED"
            except (IOError, OSError) as e:
                print(f"WARNING: log read error: {e}", file=sys.stderr)

    details = "Waiting..."
    samples_progress = []
    recent_tasks = []
    
    # If running, parse trace file
    if state == "RUNNING" and current_analysis:
        work_dir = current_analysis.get('work_dir', '')
        sample_name = current_analysis.get('sample_name', '')
        analysis_dir = current_analysis.get('analysis_dir', '')
        
        trace_file = os.path.join(analysis_dir, 'pipeline_info', 'trace.txt')
        
        if os.path.exists(trace_file):
            STAGES = {
                'ALIGN_AND_SORT': 25,
                'MARK_DUPLICATES': 30,
                'SAMTOOLS_BAM_STATS': 32,
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
                        parts = line.strip().split('	')
                        if len(parts) < 3:
                            continue
                    
                        full_name = parts[1]
                        process_name = full_name.split(' (')[0].split(':')[-1]
                        
                        pct = STAGES.get(process_name, 0)
                        if pct > progress:
                            progress = pct
            except (IOError, OSError, ValueError) as e:
                print(f"WARNING: trace parse error: {e}", file=sys.stderr)
            
            if state == 'SUCCESS':
                progress = 100
            
            samples_progress = [{
                'sample': f'{work_dir}/{sample_name}',
                'progress': progress
            }]
            
            details = f"Processing {work_dir}/{sample_name}..."
    
    return jsonify({
        'status': state,
        'details': details,
        'samples': samples_progress,
        'tasks': recent_tasks,
        'current': current_analysis
    })

@app.route('/download/<path:filepath>')
def download_result(filepath):
    """Download result file from output/."""
    # filepath format: work_dir/sample_name/subdir/filename
    target_path = os.path.join(OUTPUT_BASE, filepath)
    
    if not os.path.exists(target_path):
        return jsonify({'error': 'File not found'}), 404
    
    directory = os.path.dirname(target_path)
    filename = os.path.basename(target_path)
    
    return send_from_directory(directory, filename)

@app.route('/view_report/<work_dir>/<sample_name>')
def view_report(work_dir, sample_name):
    """Legacy URL: send users to the raw final_report.html (same as Results → view HTML). No extra wrapper page."""
    report_path = os.path.join(OUTPUT_BASE, work_dir, sample_name, 'final_report.html')
    if not os.path.exists(report_path):
        return render_template('no_report.html',
                                work_dir=work_dir,
                                sample_name=sample_name), 404
    return redirect(f'/download/{work_dir}/{sample_name}/final_report.html')

@app.route('/api/report/<work_dir>/<sample_name>')
def get_report_data(work_dir, sample_name):
    """Get report data and available files."""
    output_dir = os.path.join(OUTPUT_BASE, work_dir, sample_name)
    
    if not os.path.exists(output_dir):
        return jsonify({'error': 'Output directory not found'}), 404
    
    # List all files in output directory
    files = []
    for root, dirs, filenames in os.walk(output_dir):
        for filename in filenames:
            rel_path = os.path.relpath(os.path.join(root, filename), output_dir)
            file_size = os.path.getsize(os.path.join(root, filename))
            files.append({
                'name': filename,
                'path': rel_path,
                'size': file_size,
                'size_human': _human_readable_size(file_size)
            })
    
    return jsonify({
        'work_dir': work_dir,
        'sample_name': sample_name,
        'has_report': os.path.exists(os.path.join(output_dir, 'final_report.html')),
        'files': files
    })

def _human_readable_size(size):
    """Convert bytes to human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} TB"

@app.route('/api/logs/list', methods=['GET'])
def list_logs():
    """List all available sample logs."""
    samples = []
    
    # Scan analysis directory for *.analysis.log files
    if os.path.exists(ANALYSIS_BASE):
        for work_dir in sorted(os.listdir(ANALYSIS_BASE)):
            work_path = os.path.join(ANALYSIS_BASE, work_dir)
            if not os.path.isdir(work_path):
                continue
            
            for sample_name in sorted(os.listdir(work_path)):
                sample_path = os.path.join(work_path, sample_name)
                if not os.path.isdir(sample_path):
                    continue
                
                # Check for analysis.log file
                log_file = os.path.join(sample_path, f'{sample_name}.analysis.log')
                if os.path.exists(log_file):
                    # Check if currently running
                    is_running = (current_analysis.get('work_dir') == work_dir and 
                                current_analysis.get('sample_name') == sample_name and
                                get_pipeline_status() == "RUNNING")
                    
                    samples.append({
                        'work_dir': work_dir,
                        'sample_name': sample_name,
                        'log_file': log_file,
                        'size': os.path.getsize(log_file),
                        'modified': os.path.getmtime(log_file),
                        'is_running': is_running
                    })
    
    # Sort by modified time (most recent first)
    samples.sort(key=lambda x: x['modified'], reverse=True)
    
    return jsonify({
        'samples': samples,
        'current_sample': current_analysis if current_analysis else None
    })

@app.route('/api/logs/<work_dir>/<sample_name>', methods=['GET'])
def get_sample_log(work_dir, sample_name):
    """Get log content for a specific sample."""
    log_file = os.path.join(ANALYSIS_BASE, work_dir, sample_name, f'{sample_name}.analysis.log')
    
    if not os.path.exists(log_file):
        # Try nextflow.log as fallback
        log_file = os.path.join(LOG_BASE, work_dir, sample_name, 'nextflow.log')
        if not os.path.exists(log_file):
            return jsonify({'error': 'Log file not found'}), 404
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        return jsonify({
            'work_dir': work_dir,
            'sample_name': sample_name,
            'content': content,
            'size': os.path.getsize(log_file)
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/health')
def health():
    """Health check endpoint."""
    return jsonify({'status': 'healthy', 'service': 'dashboard'})

@app.route('/api/results')
def api_results():
    """API endpoint to list all results grouped by sample."""
    results_by_sample = {}
    
    if os.path.exists(OUTPUT_BASE):
        for work_dir in sorted(os.listdir(OUTPUT_BASE)):
            work_path = os.path.join(OUTPUT_BASE, work_dir)
            if not os.path.isdir(work_path):
                continue
            
            for sample_name in sorted(os.listdir(work_path)):
                sample_path = os.path.join(work_path, sample_name)
                if not os.path.isdir(sample_path):
                    continue
                
                sample_key = f'{work_dir}/{sample_name}'
                results_by_sample[sample_key] = {
                    'categories': {}
                }
                
                # Scan all subdirectories
                for category in sorted(os.listdir(sample_path)):
                    category_path = os.path.join(sample_path, category)
                    if not os.path.isdir(category_path):
                        continue
                    
                    files_in_category = []
                    
                    # Walk through subdirectories
                    for root, dirs, files in os.walk(category_path):
                        for filename in sorted(files):
                            filepath = os.path.join(root, filename)
                            try:
                                filesize = os.path.getsize(filepath)
                            except OSError:
                                filesize = 0
                            
                            rel_path = os.path.relpath(filepath, sample_path)
                            
                            # Determine file type and viewability
                            ext = filename.lower()
                            if ext.endswith('.html'):
                                file_type = 'HTML'
                                viewable = True
                            elif ext.endswith('.svg'):
                                file_type = 'SVG'
                                viewable = True
                            elif ext.endswith('.txt'):
                                file_type = 'Text'
                                viewable = True
                            elif ext.endswith('.vcf.gz'):
                                file_type = 'VCF.GZ'
                                viewable = False
                            elif ext.endswith('.vcf'):
                                file_type = 'VCF'
                                viewable = False
                            elif ext.endswith('.bam'):
                                file_type = 'BAM'
                                viewable = False
                            elif ext.endswith(('.bai', '.tbi')):
                                file_type = 'Index'
                                viewable = False
                            else:
                                file_type = 'Other'
                                viewable = False
                            
                            files_in_category.append({
                                'name': filename,
                                'type': file_type,
                                'path': f'{work_dir}/{sample_name}/{rel_path}',
                                'size': filesize,
                                'size_mb': round(filesize / 1024 / 1024, 2),
                                'viewable': viewable
                            })
                    
                    if files_in_category:
                        results_by_sample[sample_key]['categories'][category] = files_in_category
    
    return jsonify(results_by_sample)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
