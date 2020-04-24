# COVID19
README just for scribbling down important commands and usage

## workstation update script
`sudo echo "alias updatecovid19='curl https://raw.githubusercontent.com/KasperSkytte/covid19/workstations/singularity/updatecovid19.sh | sudo bash'" >> /etc/profile`
Then simply executing `updatecovid19` from a terminal (after a reload) will run the script with privileges
