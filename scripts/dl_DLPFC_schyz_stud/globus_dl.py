import os
from dotenv import load_dotenv
from globus_sdk import AuthClient, AccessTokenAuthorizer, NativeAppAuthClient, TransferClient, TransferData
from globus_sdk.scopes import TransferScopes

load_dotenv(".globus_env")

TRANSFER_ACCESS_TOKEN = os.getenv("ACCESS_TOKEN_TRANSFER")
TRANSFER_REFRESH_TOKEN = os.getenv("REFRESH_TOKEN_TRANSFER")
AUTH_ACCESS_TOKEN = os.getenv("ACCESS_TOKEN_AUTH")

if not TRANSFER_ACCESS_TOKEN or not TRANSFER_REFRESH_TOKEN or not AUTH_ACCESS_TOKEN:
    raise ValueError("Tokens not found in .globus_env file. Please ensure ACCESS_TOKEN and REFRESH_TOKEN are set.")
else:
    print("Tokens found, proceeding ✅")

# Create an authorizer with the access token
transfer_authorizer = AccessTokenAuthorizer(TRANSFER_ACCESS_TOKEN)
tc = TransferClient(authorizer=transfer_authorizer)

auth_authorizer = AccessTokenAuthorizer(AUTH_ACCESS_TOKEN)
ac = AuthClient(authorizer=auth_authorizer)


user = ac.oauth2_userinfo()
print("Logged in as:", user["preferred_username"])

for ep in tc.endpoint_search("Globus"):  # or any keyword you want
    print(ep["display_name"], ep["id"])

ep_id = "daf68420-f9de-42dd-a81b-3840fa43149d"

# Specify the directory to browse (you can start with root or a specific folder)
directory = "/"

items = tc.operation_ls(ep_id,  path = directory)
for item in items:
    print(item["name"], item["type"])

tc.endpoint_autoactivate(ep_id)

tdata = TransferData(
    tc,
    source_endpoint=ep_id,
    destination_endpoint="YOUR_DEST_ENDPOINT_ID",
    label="Download all FASTQs",
    sync_level="checksum"
)
for ep in tc.endpoint_search("Luxembourg"):
    print(ep["display_name"], ep["id"])