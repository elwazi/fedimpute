# Federated Imputation Workflow elwazi/fedimpute

## setting up a dev environment

1. install prerequisites:
 - nextflow
 - git-lfs
 - docker

2. clone the repo and do pip installs
```bash
# clone 
git clone git@github.com:elwazi/fedimpute.git
cd fedimpute

# create a virtual environment, and activate 
# for example direnv
echo layout python > .envrc
direnv allow .

# install the native prerequisites
# TODO: move these to a container
python -m pip install -r requirements

```

3. run the workflow with no parameters
```bash
./main.nf
```

4. or with your own vcfs
```bash 
./main.nf --input "<VCFTARGET_2.vcf.gz> <VCFTARGET_2.vcf.gz> ..etc" 
```

5. optionally, browse the output folder by serving as html
```
python -m http.server -d ./output/
```