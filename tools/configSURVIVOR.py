#!/usr/bin/env python

import argparse as ap

##### Create a custom parameter Template for Survivor ####




###########
parser = ap.ArgumentParser(description='Provide SURVIVOR runtime parameters.')
parser.add_argument('-rm','--RUNMODE',dest='run_mode',default='inv')
parser.add_argument('-svpa','--SURVPARAM', dest='sv_param', default='5000/500/5000')
###########
args = parser.parse_args()
run_mode = args.run_mode
sv_param = args.sv_param
#################################
max, min, num = map(int, sv_param.split('/'))
dup_min, dup_max, dup_num = [0, 0, 0]
ind_min, ind_max, ind_num = [0, 0, 0]
tra_min, tra_max, tra_num = [0, 0, 0]
inv_min, inv_max, inv_num = [0, 0, 0]
inv_del_min, inv_del_max, inv_del_num = [0, 0, 0]
inv_dup_min, inv_dup_max, inv_dup_num = [0, 0, 0]
####### Select RUNMODE ####
if(str.lower(run_mode) == 'dup'):
    dup_min, dup_max, dup_num = [max, min, num]
elif (str.lower(run_mode) == 'ind'):
    ind_min, ind_max, ind_num = [max, min, num]
elif (str.lower(run_mode) == 'tra'):
    tra_min, tra_max, tra_num = [max, min, num]
elif (str.lower(run_mode) == 'inv'):
    inv_min, inv_max, inv_num = [max, min, num]
############ Write template #####################
filename = 'sv_template.txt'
with open(filename, 'a') as sv_temp:
    sv_temp.write(
        'PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!\n')
    sv_temp.write('DUPLICATION_minimum_length: %d\n' % (dup_min))
    sv_temp.write('DUPLICATION_maximum_length: %d\n' % (dup_max))
    sv_temp.write('DUPLICATION_number: %d\n' % (dup_num))
    sv_temp.write('INDEL_minimum_length: %d\n' % (ind_min))
    sv_temp.write('INDEL_maximum_length: %d\n' % (ind_max))
    sv_temp.write('INDEL_number: %d\n' % (ind_num))
    sv_temp.write('TRANSLOCATION_minimum_length: %d\n' % (tra_min))
    sv_temp.write('TRANSLOCATION_maximum_length: %d\n' % (tra_max))
    sv_temp.write('TRANSLOCATION_number: %d\n' % (tra_num))
    sv_temp.write('INVERSION_minimum_length: %d\n' % (inv_min))
    sv_temp.write('INVERSION_maximum_length: %d\n' % (inv_max))
    sv_temp.write('INVERSION_number: %d\n' % (inv_num))
    sv_temp.write('INV_del_minimum_length: %d\n' % (inv_del_min))
    sv_temp.write('INV_del_maximum_length: %d\n' % (inv_del_max))
    sv_temp.write('INV_del_number: %d\n' % (inv_del_num))
    sv_temp.write('INV_dup_minimum_length: %d\n' % (inv_dup_min))
    sv_temp.write('INV_dup_maximum_length: %d\n' % (inv_dup_max))
    sv_temp.write('INV_dup_number: %d\n' % (inv_dup_num))



