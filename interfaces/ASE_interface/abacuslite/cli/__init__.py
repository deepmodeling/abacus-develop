'''abacuslite supports some complicated tasks run in cli mode
'''
import argparse
from abacuslite.cli.neb import (
    DESCRIPTION as neb_description,
    add_argument as add_neb_argument,
    run as neb_cli
)
# from abacuslite.cli.gcmc import (
#     DESCRIPTION as gcmc_description,
#     add_argument as add_gcmc_argument,
#     run as gcmc_cli
# )
DESCRIPTION = '''AbacusLite CLI - user-friendly interface for complicated tasks
Currently the cli supports following tasks:
- NEB calculation
'''

def main():
    '''
    build parser and re-direct to the corresponding implementation
    '''
    root = argparse.ArgumentParser(description=DESCRIPTION)
    subparsers = root.add_subparsers(dest='task', help='task help')

    # NEB
    neb = subparsers.add_parser('neb', help=neb_description)
    # add the argument for neb calculation
    add_neb_argument(neb)

    # # GCMC
    # gcmc = subparsers.add_parser('gcmc', help=gcmc_description)
    # # add the argument for gcmc calculation
    # add_gcmc_argument(gcmc)

    # =======================================================
    # parse
    args = root.parse_args()
    print(f'Abacuslite-CLI: task={args.task}')
    task = {
        'neb': neb_cli,
        # 'gcmc': gcmc_cli,
    }
    # run the task
    task[args.task](args)

if __name__ == '__main__':
    main()