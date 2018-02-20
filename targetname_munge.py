

def targetname_munge(target_name_in): 

    target_name_out = target_name_in 

    if 'PG-' in target_name_in[0:3]: target_name_out = 'PG '+target_name_in[3:]

    if 'WD-' in target_name_in[0:3]: target_name_out = 'WD '+target_name_in[3:]

    if 'GD-' in target_name_in[0:3]: target_name_out = 'GD '+target_name_in[3:]

    return target_name_out



 
