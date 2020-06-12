import csv
import os
import argparse
from datetime import date

FIELD_TESTS = {'SampleDate': ['date'], 'sequenced': ['yes'] }


def check_date(datestring):
    try:
        enddate_arr = datestring.split('-')  # e.g. '2020-03-22'.split('-')
        if len(enddate_arr) == 3 and len(enddate_arr[0]) == 4:
            enddate = date(year=int(enddate_arr[0]), month=int(enddate_arr[1]), day=int(enddate_arr[2]))
        elif len(enddate_arr) == 3 and len(enddate_arr[2]) == 4:
            enddate = date(year=int(enddate_arr[2]), month=int(enddate_arr[1]), day=int(enddate_arr[0]))
        else:
            enddate = None
    except ValueError as err:
        enddate = err

    return enddate


def check_file(filepath, create=False):
    """
    Checks if the supplied path contains a file, if second argument is true, creates the file if missing
    :type filepath: str
    :param filepath: relative path to check
    :type create: bool
    :param create: defaults to false , if true will create the file
    :return: absolute filepath if it does, print error message and exit otherwise
    """
    if len(filepath) == 0:
        print("invalid filepath")
        exit(-1)

    if not os.path.exists(filepath):
        if create:
            with open(filepath, 'w'):
                pass
        else:
            print("Invalid filepath: {}".format(filepath))
            exit(-1)

    return filepath


def log_field_error(field_name, row_num, err_msg, logfilewriter):
    """
    Utility to log parameter errors
    :param field_name: field where error was found
    :param row_num: row number where error was found
    :param err_msg
    :param logfilewriter: for output
    :return: None
    """
    logfilewriter.writerow(
        {'MessageType': 'Error', 'Row': row_num, 'ErrorType': 'Value Error',
         'Error in {}, details': '{}'.format(field_name, err_msg)})


def check_errors(infile, outfile, errfilewriter):
    """
    Checks for value errors, duplicates and malformed rows in metadata files
    :param infile: metadata file to check
    :param outfile: same file with value erros ommitted and
    :param errfilewriter: csv DictWriter object to log file to contain list of errors and ommitted values/rows
    :return: None
    """
    rows_read = 0
    rows_omitted = 0
    primary_keys = set()
    validated_rows = []
    with open(infile) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            rows_read += 1
            outrow = {}
            # Integrity checks
            # duplicates
            pk = row['ssi_id']
            if pk not in primary_keys:
                primary_keys.add(pk)
                outrow['ssi_id'] = pk
            else:
                err_msg = 'key: {}'.format(pk)
                errfilewriter.writerow({'MessageType': 'Error', 'Row': rows_read, 'ErrorType': 'Duplicate Key',
                                        'Details': '{}'.format(err_msg)})
                rows_omitted += 1
                continue

            # Field checks
            for field_name in row.keys():
                if field_name in FIELD_TESTS:
                    val: str = row[field_name]
                    for test in FIELD_TESTS[field_name]:
                        if test == 'date':
                            res = check_date(val)
                            if isinstance(res, ValueError):
                                log_field_error(field_name, rows_read, "Invalid date value: {}".format(val),
                                                errfilewriter)
                                outrow[field_name] = ''
                            elif isinstance(res, date):
                                outrow[field_name] = res.strftime("%Y-%m-%d")
                        if test == 'yes':
                            if len(val)>0:
                                if val.lower() in ['y','yes','ja','true','sand']:
                                    outrow[field_name] = 1
                                elif val.lower() in ['n','no','nej','false','falsk']:
                                    outrow[field_name] = 0
                                else:
                                    log_field_error(field_name, rows_read, "Invalid yes/no value: {}".format(val)
                                                    , errfilewriter)

            # TODO number of fields in row
            # TODO Common sense checks

            validated_rows.append(outrow)

    print("Finished processing {} rows of which {} where ommitted due to irrecoverable errors".format(rows_read,
                                                                                                      rows_omitted))
    with open(outfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, reader.fieldnames)
        writer.writeheader()
        for row in validated_rows:
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='review a metadata file to check for errors')
    parser.add_argument('infile', type=str, help='path to the file to be checked')
    parser.add_argument('outfile', type=str, help='path to the clean file to be created')
    parser.add_argument('errfile', type=str, help='path to the error file to be created')
    args = parser.parse_args()
    infile = check_file(args.infile)
    print('Validated infile as {}'.format(infile))
    outfile = check_file(args.outfile, True)
    print('Validated outfile as {}'.format(outfile))
    errfile = check_file(args.errfile, True)
    print('Validated errfile as {}'.format(errfile))
    with open(errfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, ['MessageType', 'Row', 'ErrorType', 'Details'])
        writer.writeheader()
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Started {}'.format(infile)})
        check_errors(infile, outfile, writer)
        writer.writerow({'MessageType': 'Info', 'ErrorType': '', 'Details': 'Finished {}'.format(infile)})
