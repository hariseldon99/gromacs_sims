
## Antimicrobial_Peptides
- Contains subdirectories for various peptide simulations.
- Anti-Microbial Research Association (AMRA), The University of Burdwan 
  - [Professor Rajib Bandopadhyay](http://rbandopadhyayslab.unaux.com/), Department of Botany 
  - [Dr. Analabha Roy](https://physics.utexas.edu/~daneel), Department of Physics
  - [Dr. Sumit Hira](https://www.sumithira.in/), Department of Zoology
  - Rajendra Kr Roy, Department of Botany 
  - Dr. Raju Biswas, Department of Botany  
  - Rajdeep Shaw, Department of Botany  
  - Dr. Rahul Chandra, Department of Physics 
  - Samrat Daripa, Department of Zoology
  - Argha Nath, Department of Zoology

- Contents:
  - `Colpk_phosphatedylethanolamine/`, `Colpk_phosphatedylglycerol/`, `Colpk_pyocyanin/`: Simulations involving Colicin peptides.
  - `Colpm_phosphatedylethanolamine/`, `Colpm_phosphatedylglycerol/`, `Colpm_pyocyanin/`: Simulations involving Colicin PM peptides.
  - `DNA_peptide/`: Simulations involving DNA-peptide interactions.
  - `Ku04AMP01_linear/`, `Ku04AMP01_phosphatedylethanolamine/`: Simulations of Ku04AMP01 antimicrobial peptide.

- Deletion Note on 20250505
There was a google drive synchroniztion error that led to topology and trajectory data loss. Full restoration was accomplished on this day. Restoration script:
```python
# pip install --upgrade google-api-python-client google-auth-httplib2 google-auth-oauthlib
from __future__ import print_function
import os
import pickle
import datetime
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

# If modifying these SCOPES, delete the file token.pickle.
SCOPES = ['https://www.googleapis.com/auth/drive']

def authenticate():
    """
    Handles user authentication. If valid credentials are stored in token.pickle, they are used;
    otherwise, the OAuth flow starts using the downloaded credentials.json file.
    """
    creds = None
    # Token file stores the user's credentials.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            # credentials.json should be your downloaded OAuth 2.0 client credentials.
            flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run.
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
    return creds

def main():
    # Authenticate and build the Drive API service.
    creds = authenticate()
    service = build('drive', 'v3', credentials=creds)

    # Compute the RFC 3339 timestamp for three months (90 days) ago from now (UTC).
    four_months_ago = datetime.datetime.utcnow() - datetime.timedelta(days=120)
    four_months_ago_str = four_months_ago.isoformat() + 'Z'

    # Build a query to find files in the trash that:
    # - have "trashed" set to true,
    # - whose name contains ".xtc" or ".tpr", and
    # - whose modifiedTime is later than three months ago.
    #
    # Note: Google Drive does not provide a dedicated "deletion date" field,
    # so we use modifiedTime as a proxy. Depending on your usage, this may not perfectly reflect
    # when the file was trashed.
    query = (
        f"trashed = true and "
        #f"(name contains '.xtc' or name contains '.tpr') and "
        f"(name contains '.xtc' or name contains 'colpk') and "
        f"modifiedTime > '{four_months_ago_str}'"
    )

    files_to_restore = []
    page_token = None
    print("Searching for trashed files to restore...\n")
    while True:
        response = service.files().list(
            q=query,
            spaces='drive',
            fields='nextPageToken, files(id, name)',
            pageToken=page_token
        ).execute()
        for file in response.get('files', []):
            print(f"Found file: {file.get('name')} (ID: {file.get('id')})")
            files_to_restore.append(file)
        page_token = response.get('nextPageToken', None)
        if page_token is None:
            break

    if not files_to_restore:
        print("\nNo matching files found in the trash.")
        return

    # Restore each file by updating its 'trashed' property to False.
    print("\nRestoring files...")
    for file in files_to_restore:
        file_id = file.get('id')
        updated_file = service.files().update(
            fileId=file_id,
            body={'trashed': False}
        ).execute()
        print(f"Restored file: {updated_file.get('name')} (ID: {updated_file.get('id')})")

if __name__ == '__main__':
    main()
```
**Note:** Get google auth credentials.json by the following steps:

Go to the Google Cloud Console: Open your browser and navigate to the Google Cloud Console. Make sure you’re logged in with the account that owns the project where you enabled the Drive API.

Select Your Project: From the top navigation, select the project that has the Google Drive API enabled.

Open the Credentials Page: In the left sidebar, click on "APIs & Services", then choose "Credentials."

Configure the OAuth Consent Screen (if not done already): If you haven’t already set this up, you’ll see a prompt or an option to configure the OAuth consent screen. Click on it and follow the instructions. You’ll need to provide some basic details like the application name, support email, and optionally scopes and authorized domains. This configuration is necessary before creating an OAuth client.

Create OAuth Client Credentials: Click on "Create Credentials" at the top of the Credentials page and select "OAuth client ID".

For Desktop Applications: If you are running the script locally, choose the “Desktop app” option.

For Web Applications: If you plan to run it on a server, choose “Web application” and fill in the required fields (like authorized redirect URIs).

Download the Credentials File: Once you create the OAuth client ID, a dialog box will display your client details. Look for the "Download JSON" button. Clicking it downloads the file (commonly named credentials.json) onto your computer.

Place the File in Your Project: Save the downloaded credentials.json in the working directory of your Python project so that your script can locate it during authentication.