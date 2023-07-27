import os,json,traceback,sys

# This script will read the machine.json file and modify the username, password and project id

username = os.environ.get("BOHRIUM_USERNAME",None)
password = os.environ.get("BOHRIUM_PASSWORD",None)
projectid = os.environ.get("BOHRIUM_PROJECT_ID",None)

try:
    machine_file = json.load(open("machine.json"))
    if username:
        machine_file["train"][0]["machine"]["remote_profile"]["email"] = username
        machine_file["model_devi"][0]["machine"]["remote_profile"]["email"] = username
        machine_file["fp"][0]["machine"]["remote_profile"]["email"] = username

    if password:
        machine_file["train"][0]["machine"]["remote_profile"]["password"] = password
        machine_file["model_devi"][0]["machine"]["remote_profile"]["password"] = password
        machine_file["fp"][0]["machine"]["remote_profile"]["password"] = password

    if projectid:
        projectid = int(projectid)
        machine_file["train"][0]["machine"]["remote_profile"]["program_id"] = projectid
        machine_file["model_devi"][0]["machine"]["remote_profile"]["program_id"] = projectid
        machine_file["fp"][0]["machine"]["remote_profile"]["program_id"] = projectid

    if os.environ.get("ABACUS_IMAGE"):
        machine_file["fp"][0]["machine"]["remote_profile"]["input_data"]["image_name"] = os.environ.get("ABACUS_IMAGE")

    json.dump(machine_file,open("machine.json","w"),indent=4)
except:
    traceback.print_exc()
    sys.exit(1)

sys.exit(0)
