#!/usr/bin/env python
import io
import argparse
import requests
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build

# Google Drive API setup
SCOPES = ['https://www.googleapis.com/auth/drive.file']

def authenticate_google_drive(credentials_path):
    creds = None
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)
    else:
        flow = InstalledAppFlow.from_client_secrets_file(credentials_path, SCOPES)
        creds = flow.run_local_server(port=0)
        with open('token.json', 'w') as token:
            token.write(creds.to_json())
    return build('drive', 'v3', credentials=creds)

# Function to stream a file from Google Drive
def stream_file(service, file_id):
    request = service.files().get_media(fileId=file_id)
    return request

# Zenodo upload
def upload_to_zenodo(file_stream, file_name, access_token):
    headers = {"Authorization": f"Bearer {access_token}"}
    data = {
        'metadata': {
            'title': 'My Upload',
            'upload_type': 'dataset'
        }
    }
    response = requests.post('https://zenodo.org/api/deposit/depositions', headers=headers, json=data)
    deposition_id = response.json()['id']
    bucket_url = response.json()['links']['bucket']

    response = requests.put(f"{bucket_url}/{file_name}", data=file_stream, headers=headers)
    if response.status_code == 200:
        print("Upload successful!")
    else:
        print("Upload failed!")

# Main function
def main():
    parser = argparse.ArgumentParser(description='Stream a file from Google Drive and upload it to Zenodo.')
    parser.add_argument('google_credentials', help='Path to the Google credentials JSON file')
    parser.add_argument('zenodo_api_key', help='Zenodo API key')
    parser.add_argument('file_id', help='Google Drive file ID to stream. Right-click the file name and select "Get shareable link". The last part of the link is the file ID')
    parser.add_argument('file_name', help='Name of the file to upload to Zenodo')

    args = parser.parse_args()

    drive_service = authenticate_google_drive(args.google_credentials)
    file_stream = stream_file(drive_service, args.file_id)
    upload_to_zenodo(file_stream, args.file_name, args.zenodo_api_key)

if __name__ == '__main__':
    main()