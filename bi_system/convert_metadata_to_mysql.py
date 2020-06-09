import mysql.connector

host = input("Please enter the mysql database server IP (defaults to localhost): ")
host = 'localhost' if len(host.strip(' ')) == 0 else host
password = input("Please enter the mysql user's database password: ")
cnx = mysql.connector.connect(user='covid19', password=password,
                              host=host,
                              database='covid19')
print('success')
cnx.close()
