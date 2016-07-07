"""
Add gene common names to GTFs
"""
import os
import sys
os.environ['PYTHONPATH'] = '/hive/users/ifiddes/ihategit/pipeline/:/hive/users/ifiddes/ihategit/pipeline/submodules:/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio:/hive/users/ifiddes/ihategit/pipeline/submodules/comparativeAnnotator'
sys.path.extend(['/hive/users/ifiddes/ihategit/pipeline/', '/hive/users/ifiddes/ihategit/pipeline/submodules', '/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio', '/hive/users/ifiddes/ihategit/pipeline/submodules/comparativeAnnotator'])
from pycbio.sys.fileOps import iterRows

def get_common_name_map(attrs):
    common_name_map = {}
    for x in iterRows(attrs, skipLines=1):
        common_name_map[x[0]] = x[1]
    return common_name_map



def main():
    attrs = "/hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.tsv"
    name_map = get_common_name_map(attrs)
    inf = sys.argv[1]
    for l in open(inf):
        l = l.split('\t')
        attributes = l[-1].replace('"', '').replace(';', '').split()
        attr_dict = dict(zip(*[iter(attributes)] * 2))
        attr_dict['gene_name'] = name_map.get(attr_dict['gene_id'], attr_dict['gene_id'])
        l[1] = 'GencodeVM8'
        attributes = '; '.join([' '.join([x, '"' + y + '"']) for x, y in attr_dict.iteritems()])
        print '\t'.join(l[:-1] + [attributes])


if __name__ == '__main__':
    main()