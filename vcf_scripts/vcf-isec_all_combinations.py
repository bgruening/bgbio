#/usr/bin/env python

import os
import sys
import itertools
import subprocess
import argparse
import tempfile

def simple_line_count( filename ):
    """
        Can be easily optimized.
    """
    lines = 0
    for not line in open(filename):
        if line.startswith('#'):
            lines += 1
    return lines

def get_all_combinations( path_list ):
    """
        Calculates all possible combinations of a given list
    """
    for combi_len in range(2, len(path_list)+1):
        for combinations in itertools.combinations( path_list, combi_len ):
            yield combinations

def run_vcf_intersect( paths ):
    """
        runs vcf_isec on a list of given paths
    """
    command = 'vcf-isec -f %s' % ( ' '.join( paths ) )
    #use for simple testing
    #command = 'cat %s' % ( ' '.join( paths ) )
    stderr_name = tempfile.NamedTemporaryFile( prefix = "vcf_isec_stderr" ).name
    stdout_name = tempfile.NamedTemporaryFile( prefix = "vcf_isec_stdout" ).name
    try:
        subprocess.check_call( args=command, shell=True, stderr=open( stderr_name, 'wb' ), stdout=open( stdout_name, 'wb' ) )
        return simple_line_count( stdout_name )
    except:
        for line in open( stderr_name ):
            sys.stderr.write( line )
        raise Exception( "Error intersecting VCF files" )
    finally:
        os.unlink( stderr_name )


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Runs vcf-isec on every combination of given vcf files.')

    parser.add_argument("-i", "--inputs", nargs="*",
                    required=True,
                    help="Path to the VCF input files.")

    args = parser.parse_args()
    if len(args.inputs) < 2:
        sys.exit('Minimum 2 files needs to be given.')
    combinations = get_all_combinations( args.inputs )
    for combination in combinations:
        common_lines = run_vcf_intersect( combination )
        print '%s\t%s' % ('\t'.join(combination), common_lines)

