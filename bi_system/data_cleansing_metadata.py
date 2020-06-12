import csv
import os
import argparse
from datetime import date

FIELD_TESTS = dict(SampleDate=['date'], sequenced=['yes'], Sex=['vals:F_M'], Travel=['yes'], ContactWithCase=['yes'],
                   EpilprResp_start=[date], EpilprResp=['yes'], EpilprVent_start=[date], EpilprVent=['yes'],
                   EpilprECMO_start=[date], EpilprECMO=['yes'], EpilprHeart=['yes'], EpilprHeart_start=[date],
                   Diabet=['yes'], Neuro=['yes'], Cancer=['yes'], Adipos=['yes'], Nyre=['yes'], Haem_c=['yes'],
                   Card_dis=['yes'], Resp_dis=['yes'], Immu_dis=['yes'], Other_risk=['yes'], Pregnancy=['yes'],
                   Doctor=['yes'], Nurse=['yes'], PlaceOfInfection_EN=['dim:countries.tsv'], ReportAge=['age'],
                   ReportAgeGrp=['dim:age_groups.tsv'], COVID19_Status=['vals:0_1_2:0'], COVID19_EndDate=['date'],
                   lineage=['str'], lineages_version=['date'], ParishCode=['dim:parish.tsv'],
                   MunicipalityCode=['dim:municipalities.tsv'], NUTS3Code=['dim:nuts3_regions.csv'], Occupation=['str'],
                   CountryOfTravel=['dim:countries.tsv'], SymptomsStartDate=['date'], CodR_DateOfDeath=['date'],
                   CodR_Death60Days=['yes'], CPR_Death60Days=['yes'], CPR_DateOfDeath=['date'],
                   DateOfDeath_final=['date'])

"""
date: check that date is well-formed and output in 'dddd-mm-yy' format
yes: check that field is y/n or 1/0 or any other boolean variation and output as 1/0
vals:val1_val2_...[:default]: Check that val is in list of given values. 
                              Optionally output default value when no value is supplied
dim:dim_name: check that code exists in first column of stable dimension file.
age: check int and reasonable number (0-120)
str: string - just pass through 
"""


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
         'Details': 'Error in {}, details: {}'.format(field_name, err_msg)})


def load_dims():
    """
    Load keys of dimension from dimension file
    :return: dictionary from dimension filename to set of keys
    """
    dims = dict()
    for test_list in FIELD_TESTS.values():
        for test in test_list:
            if test.startswith('dim'):
                dim_filename = test.split(':')[1]
                dim_name = dim_filename.split('.')[0]
                if not dim_name in dims.keys():
                    dim_keys = set()
                    with open(os.path.join("stable_dims", dim_filename)) as csvfile:
                        reader = csv.reader(csvfile, delimiter='\t' if dim_filename.endswith('tsv') else ',')
                        next(reader, None)  # skip header
                        for row in reader:
                            if len(row[0]) > 0:
                                dim_keys.add(row[0])
                    dims[dim_name] = dim_keys

    return dims


def check_errors(datafile, outfile, errfilewriter):
    """
    Checks for value errors, duplicates and malformed rows in metadata files
    :param datafile: metadata file to check
    :param outfile: same file with value erros ommitted and
    :param errfilewriter: csv DictWriter object to log file to contain list of errors and ommitted values/rows
    :return: None
    """

    static_dims = load_dims()
    rows_read = 0
    rows_omitted = 0
    primary_keys = set()
    validated_rows = []
    with open(datafile) as csvfile:
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
                            if len(val) > 0:
                                if val.lower() in ['y', 'yes', 'ja', 'true', 'sand', '1']:
                                    outrow[field_name] = 1
                                elif val.lower() in ['n', 'no', 'nej', 'false', 'falsk', '0']:
                                    outrow[field_name] = 0
                                else:
                                    log_field_error(field_name, rows_read, "Invalid yes/no value: {}".format(val)
                                                    , errfilewriter)
                        if test.startswith('vals'):
                            test_params = test.split(':')
                            allowed_vals = test_params[1].split('_')
                            if len(val) > 0:
                                if val in allowed_vals:
                                    outrow[field_name] = val
                                else:
                                    log_field_error(field_name, rows_read, "Invalid value: {}, expected one of {}"
                                                    .format(val, allowed_vals), errfilewriter)
                            else:
                                if len(test_params) == 3:  # there exists a default value
                                    outrow[field_name] = test_params[2]

                        if test.startswith('dim'):
                            dim_name = test.split(':')[1].split('.')[0]
                            keys = static_dims[dim_name]
                            if len(val) > 0:
                                if not val in keys:
                                    if val == 'Not Denmark, Unknown':
                                        outrow[field_name] = val  # pass through
                                    log_field_error(field_name, rows_read, "Invalid value: {}, expected corresponding"
                                                                           "dimension key in {}"
                                                    .format(val, dim_name), errfilewriter)
                                else:
                                    outrow[field_name] = val

                        if test == 'age':
                            if len(val) > 0:
                                try:
                                    int_val = int(val)
                                    if int_val >= 0 and int_val < 120:
                                        outrow[field_name] = val
                                    else:
                                        log_field_error(field_name, rows_read,
                                                        "Invalid numeric value: {}, expected between 0 and 120"
                                                        .format(val), errfilewriter)
                                except Exception as e:
                                    log_field_error(field_name, rows_read, "Non numeric value: {}"
                                                    .format(e), errfilewriter)

                        if test == 'str':  # just pass through
                            outrow[field_name] = val

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
