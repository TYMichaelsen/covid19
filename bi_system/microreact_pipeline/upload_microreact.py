from google.oauth2 import service_account
from googleapiclient.discovery import build
from apiclient.http import MediaFileUpload
from apiclient import errors

METADATA_DRIVE_FILE_ID = '11_5qlSGVcYGb-x8J30igO9_kXlC-UGUO'
TREE_DRIVE_FILE_ID = '17J9ke-_KVmj3vhQq_gXdF6WME_LJG1pX'


def upload(config):    
    credentials = _get_credentials()
    service = build('drive', 'v3', credentials=credentials)
    _upload_file(config['out_react_tsv'], METADATA_DRIVE_FILE_ID, service)
    _upload_file(config['out_react_nwk'], TREE_DRIVE_FILE_ID, service)

def _get_credentials():
    try:
        scopes = ['https://www.googleapis.com/auth/drive']
        service_acc_file = './microreact_pipeline/config/service_account.json'
        return service_account.Credentials.from_service_account_file(service_acc_file, scopes=scopes)
    except Exception as e:
        print(e)
        print('Failed to obtain credentials to google service account maintaining microreact file.')

def _upload_file(new_filepath, drive_file_id, service):
     try:
        new_file_content = MediaFileUpload(new_filepath, resumable=True)
        result = service.files().update(fileId=drive_file_id, media_body=new_file_content).execute()
        if result is None:
            raise Exception('Google drive could not update file with ID {}'.format(drive_file_id))
    except errors.HttpError as e:
        print(e)
        print('Error occured when attempting to update {} on google drive.'.format(drive_file_id))