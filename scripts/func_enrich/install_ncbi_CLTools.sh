# This script installs datasets and dataformat NCBI Command Line tools in a ~/.local/bin directory
# (within user's home folder). If this directory does not exist it creates it.

# Check if the directory exists, and if not create it
locBinDir=~/.local/bin
if [ ! -d $locBinDir ]
then
    mkdir $locBinDir
    cd $locBinDir
else
    cd $locBinDir
fi

# If the .local/bin directory is not in path, add it to the .bashrc file and
# load it again
if [[ ! ":$PATH:" == *":${locBinDir}:"* ]]
then
    newPath="${locBinDir}:$PATH"
    echo "${locBinDir} is not in path. Adding it."
    echo "export PATH=$newPath" >> ~/.bashrc
    source ~/.bashrc
else
    echo "${locBinDir} is already in path."
fi

# Check if the tools are already downloaded and, if not, download them
if [ ! -f datasets ]
then
    curl https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets --output datasets
    chmod 777 datasets
    echo "datasets downloaded successfully and installed at ${locBinDir}"
else
    echo "datasets already exists in ${locBinDir}"
fi

if [ ! -f dataformat ]
then
    curl https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat --output dataformat
    chmod 777 dataformat
    "dataformat downloaded successfully and installed at ${locBinDir}"
else
    echo "dataformat already exists in ${locBinDir}"
fi
