#!/usr/bin/env python3
"""Verify the archived resource and reproduce the article's numerical tables."""
import argparse,json,platform,sys,time,subprocess,shutil
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/'scripts'))
import verify_release
import summarize_resource

def argument_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',required=True,type=Path,help='Unmodified resource export directory')
    p.add_argument('--out',required=True,type=Path,help='New results directory')
    followup=p.add_mutually_exclusive_group()
    followup.add_argument('--followup',type=Path,help='Completed follow-up read-audit archive (optional)')
    followup.add_argument('--sd13',dest='followup',type=Path,help=argparse.SUPPRESS)
    masks=p.add_mutually_exclusive_group()
    masks.add_argument('--masks',type=Path,help='Archived common-observation masks (optional)')
    masks.add_argument('--maskdata',dest='masks',type=Path,help=argparse.SUPPRESS)
    return p

def main():
    p=argument_parser()
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
        if a.followup is not None:
            import verify_followup
            # Preserve recorded output keys for existing downstream readers.
            status['sd13']=verify_followup.verify(a.followup,out/'targeted_audit')
        if a.masks is not None:
            import verify_recovered_mask
            status['recovered_mask']=verify_recovered_mask.verify_recovered_mask(a.masks,out/'recovered_mask')
        status.update(status='PASS',elapsed_seconds=round(time.monotonic()-start,3),verification=verification,summary=report)
    except Exception as e:
        status.update(status='FAILED',error=repr(e));raise
    finally:
        (out/'run_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print('RESOURCE_REPLAY=PASS')
    print('OUTPUT='+str(out))
if __name__=='__main__':main()
