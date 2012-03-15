import macroeco.params as params

pa = params.Parameters('callparams', {'size':'number', 'species':'foo','layers':'list of integers'})
print pa.params
print pa.interactive
onerun = pa.params['autoname0']
print map(int, onerun['layers'][1:-1].split(','))
