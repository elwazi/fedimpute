#!/usr/bin/env python3

import requests
import json
import argparse
import sys
import time
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--legend", help="legend file")
parser.add_argument("--hap", help="hap file")
parser.add_argument(
    "--outfile", help="output file pattern (will end with .hap.gz and .legend.gz)")
parser.add_argument("--start",  type=int, help="start position")
parser.add_argument("--stop",  type=int, help="stop position")
args = parser.parse_args()


def submit_job(input_files, hap_file, outfile, start, stop):

    # imputation server url
    url = 'https://imputationserver.sph.umich.edu/api/v2'

    # add token to header (see Authentication)
    headers = {'X-Auth-Token': token}

    # submit new job
    vcf = '/path/to/genome.vcf.gz'
    files = {'input-files': open(vcf, 'rb')}
    r = requests.post(url + "/jobs/submit/minimac4",
                      files=files, headers=headers)
    if r.status_code != 200:
        raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))

    # print message
    print r.json()['message']
    print r.json()['id']

    return None


if __name__ == '__main__':
    submit_job()
