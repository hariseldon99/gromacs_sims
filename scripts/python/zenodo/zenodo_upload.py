import os
import sys
import socket
import argparse
import urllib3
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from tqdm import tqdm

# Parse arguments
parser = argparse.ArgumentParser(description="Upload a file to a Zenodo deposition.")
parser.add_argument("-d", "--deposition-id", required=True, help="Zenodo deposition ID")
parser.add_argument("-f", "--filepath", required=True, help="Full path to the file to upload")
args = parser.parse_args()

# Configuration
ACCESS_TOKEN = os.environ.get("ZENODO_API_TOKEN")
if not ACCESS_TOKEN:
    print("Error: ZENODO_API_TOKEN environment variable is not set.")
    sys.exit(1)

DEPOSITION_ID = args.deposition_id
FILEPATH = args.filepath
FILENAME = os.path.basename(FILEPATH)


def human_readable_size(size_bytes):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"


# 1. Setup a resilient connection session
session = requests.Session()
retries = Retry(total=10, backoff_factor=2, status_forcelist=[500, 502, 503, 504])
session.mount('https://', HTTPAdapter(max_retries=retries))
session.headers.update({"Authorization": f"Bearer {ACCESS_TOKEN}"})

print("Fetching bucket URL...")
resp = session.get(f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}")
resp.raise_for_status()
bucket_url = resp.json()['links']['bucket']
target_url = f"{bucket_url}/{FILENAME}"

file_size = os.path.getsize(FILEPATH)
print(f"File      : {FILENAME}")
print(f"Full size : {human_readable_size(file_size)} ({file_size:,} bytes)")

# 2. Check for an existing partial upload and resume from that offset
offset = 0
head_resp = session.head(target_url)
if head_resp.status_code == 200:
    remote_size = int(head_resp.headers.get('Content-Length', 0))
    if remote_size >= file_size:
        print("File already fully uploaded. Nothing to do.")
        sys.exit(0)
    elif remote_size > 0:
        offset = remote_size
        print(f"Resuming  : {human_readable_size(offset)} already uploaded, continuing from byte {offset:,}")

remaining = file_size - offset
print(f"Remaining : {human_readable_size(remaining)} ({remaining:,} bytes)")
print(f"Target    : {target_url}")
print("Speed and ETA will be shown in the progress bar. Please leave your machine awake.")

# 3. Enable socket-level keep-alives to prevent Zenodo from dropping quiet links
session.get_adapter('https://').poolmanager.connection_pool_kw['socket_options'] = [
    (socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1),
]

# 4. Stream the file from the resume offset
upload_headers = {}
if offset > 0:
    upload_headers['Content-Range'] = f"bytes {offset}-{file_size - 1}/{file_size}"
    upload_headers['Content-Length'] = str(remaining)

try:
    with open(FILEPATH, 'rb') as f:
        f.seek(offset)
        with tqdm.wrapattr(
            f, 'read',
            total=remaining,
            initial=0,
            desc=FILENAME,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
            dynamic_ncols=True,
        ) as fobj:
            upload_resp = session.put(target_url, data=fobj, headers=upload_headers)
    upload_resp.raise_for_status()
    print(f"\nSuccess! Uploaded {human_readable_size(file_size)} to deposition {DEPOSITION_ID}.")
except Exception as e:
    print(f"\nError: {e}")
    sys.exit(1)