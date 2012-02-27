import params

pa = params.Parameters('Unnamed', {'size':'number', 'species':'foo','layers':'list of integers'})
print pa.params
print pa.interactive
onerun = pa.params['autoname0']
print map(int, onerun['layers'][1:-1].split(','))
