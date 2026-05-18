import argparse

from pandda_gemmi.pandda.ligand_build import main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    
    parser.add_argument('--dataset_dir')
    args = parser.parse_args()
    print(args)

    main(
        args.dataset_dir,
    )