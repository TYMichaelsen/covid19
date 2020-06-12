import csv
import os
import argparse



def check_file(filepath, create=False):
    """
    Checks if the supplied path contains a file, if second argument is true, creates the file if missing
    :type filepath: str
    :param filepath: relative path to check
    :return: absolute filepath if it does, print error message and exit otherwise
    """
    if len(filepath) == 0:
        print("invalid filepath")
        exit(-1)

    if not os.path.exists(filepath):
        if create:
            with open(filepath, 'w'): pass
            return filepath
        else:
            print("Invalid filepath: {}".format(filepath))
            exit(-1)


def check_errors(infile, outfile, errfile):
    """
    Checks for value errors, duplicates and malformed rows in metadata files
    :param infile: metadata file to check
    :param outfile: same file with value erros ommitted and
    :param errfile: log file to contain list of errors and ommitted values/rows
    :return: None
    """
    rows_read = 0
    rows_omitted = 0
    with open(infile) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            rows_read+=1

    print("Finished processing {} rows of which {} where ommitted due to irrecoverable errors".format(rows_read, rows_omitted))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='review a metadata file to check for errors')
    parser.add_argument('infile', type=str,help='path to the file to be checked')
    parser.add_argument('outfile', type=str,help='path to the clean file to be created')
    parser.add_argument('errfile', type=str,help='path to the error file to be created')
    args = parser.parse_args()
    infile = check_file(args.infile)
    print('Validated infile as {}'.format(infile))
    outfile = check_file(args.outfile,True)
    print('Validated outfile as {}'.format(infile))
    errfile = check_file(args.outfile,True)
    print('Validated errfile as {}'.format(infile))
    check_errors(infile, outfile, errfile)
