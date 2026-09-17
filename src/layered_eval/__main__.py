"""Portable CLI: prepare raw counts, validate bundle, evaluate measured responses."""
import argparse,json
from pathlib import Path
from .bundle import validate_bundle
from .pipeline import prepare_counts,evaluate_measured
from .predictions import score_predictions

def main():
 p=argparse.ArgumentParser(description='Layered perturbation evidence evaluation')
 sub=p.add_subparsers(dest='command',required=True)
 v=sub.add_parser('validate');v.add_argument('bundle');v.add_argument('--output')
 a=sub.add_parser('prepare');a.add_argument('config');a.add_argument('--output',required=True)
 e=sub.add_parser('evaluate');e.add_argument('bundle');e.add_argument('--output',required=True);e.add_argument('--iterations',type=int,default=100);e.add_argument('--n',type=int,default=16)
 s=sub.add_parser('score-predictions');s.add_argument('config');s.add_argument('--output',required=True)
 args=p.parse_args()
 if args.command=='validate':
  r=validate_bundle(args.bundle)
  if args.output:Path(args.output).write_text(json.dumps(r,indent=2)+'\n')
 elif args.command=='prepare':r=prepare_counts(args.config,args.output)
 elif args.command=='score-predictions':r=score_predictions(args.config,args.output)
 else:r=evaluate_measured(args.bundle,args.output,args.iterations,args.n)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
