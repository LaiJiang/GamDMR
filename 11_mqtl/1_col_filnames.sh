python3 - <<'PY'
import requests

url = "https://zenodo.org/api/records/8047777"

r = requests.get(url, timeout=60)
r.raise_for_status()
j = r.json()

for f in j["files"]:
    print(
        f["key"],
        round(f["size"]/1024**3, 3), "GB",
        f["links"].get("content", f["links"].get("self"))
    )
PY


#######after collecting file names, download the files with wget, e.g.:

cd /mnt/c/Per/LaiJiang/Project/UQAC/meth/dat/EPIGEN

wget -c -O README.md \
"https://zenodo.org/api/records/8047777/files/README.md/content"

wget -c -O meQTL_full.txt.gz \
"https://zenodo.org/api/records/8047777/files/meQTL_full.txt.gz/content"


#after download, verify the md5sum of the files
md5sum meQTL_full.txt.gz

gzip -t meQTL_full.txt.gz
echo $?