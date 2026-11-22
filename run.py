#!/usr/bin/env python3
"""Verify the frozen export, then reproduce manuscript source tables."""
import argparse,json,platform,sys,time,subprocess,shutil
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/'scripts'))
import verify_release
import summarize_resource

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',required=True,type=Path,help='Unmodified release_candidate directory')
    p.add_argument('--out',required=True,type=Path,help='New results directory')
    p.add_argument('--sd13',type=Path,help='Completed targeted-audit REVIEW run directory (optional)')
    p.add_argument('--maskdata',type=Path,help='Recovered frozen common-mask files (optional)')
    a=p.parse_args()
    if sys.version_info<(3,9):p.error('Python 3.9+ required')
    if sys.flags.optimize:p.error('Do not run with -O')
    if a.out.exists():p.error('Output directory already exists; choose a new path')
    data=a.data.resolve();out=a.out.resolve();out.mkdir(parents=True)
    start=time.monotonic()
    status={'status':'RUNNING','python':sys.version,'platform':platform.platform()}
    try:
        verification=verify_release.check(data,out/'verification')
        report=summarize_resource.summarize(data,out/'tables')
        subprocess.run([sys.executable,str(Path(__file__).resolve().parent/'scripts/summarize_panel_structure.py'),'--data',str(data),'--out',str(out/'panel_structure')],check=True)
        for name in ('structure_recount.json','primary142_isolate_call_counts.tsv'):
            shutil.copy2(out/'panel_structure'/name,out/'tables'/name)
        if a.sd13 is not None:
            import replay_sd13
            status['sd13']=replay_sd13.replay(a.sd13,out/'targeted_audit')
        if a.maskdata is not None:
            import verify_recovered_mask
            status['recovered_mask']=verify_recovered_mask.verify_recovered_mask(a.maskdata,out/'recovered_mask')
        status.update(status='PASS',elapsed_seconds=round(time.monotonic()-start,3),verification=verification,summary=report)
    except Exception as e:
        status.update(status='FAILED',error=repr(e));raise
    finally:
        (out/'run_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print('RESOURCE_REPLAY=PASS')
    print('OUTPUT='+str(out))
if __name__=='__main__':main()
