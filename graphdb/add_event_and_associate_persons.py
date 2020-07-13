from py2neo import Graph, Node, Relationship
import argparse
import csv
import os


def log_field_error(field_name: str, row_num: int, err_msg: str, logfilewriter: csv.DictWriter) -> None:
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


def connect():
    host = input("Please enter the neo4j database server IP (defaults to localhost): ")
    host = 'localhost' if len(host.strip(' ')) == 0 else host
    password = input("Please enter the neo4j user's database password: ")
    graph = Graph("bolt://{}:7687".format(host), user='neo4j',
                  password=password)  # this may need to change when running from the server since there it needs the instance IP
    graph.delete_all()
    return graph


def load_data(graph, datafile, eventname):
    tx = graph.begin()

    # Create event
    n = Node("Event", name=eventname)
    tx.create(n)

    # Create relationships with persons
    with open(datafile) as file:
        reader = csv.reader(file)
        i = 0
        for row in reader:
            i += 1
            if len(row) > 0 and len(row[0]) > 0:
                person = graph.find_one(label='Person', property_key='ssi_id', property_value=row[0])
                if person is not None:
                    tx.create(Relationship(person, "ASSOCIATED_WITH", n))

        tx.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='load a text file with a list of cases to neo4j')
    parser.add_argument('infile', type=str, help='path to the list of SSI-IDs in format SSI-1234 and no header row')
    parser.add_argument('event', type=str, help='name of the event to be created')
    args = parser.parse_args()
    infile = check_file(args.infile)
    print('Validated infile as {}'.format(infile))
    print('Creating event {}'.format(args.event))
    graph = connect()
    load_data(graph, infile, args.event)
