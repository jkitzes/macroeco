# -*- mode: python -*-
a = Analysis(['desktop.py'],
             pathex=['/Users/jkitzes/Projects/macroeco/macroeco'],
             hiddenimports=['scipy.special._ufuncs_cxx'],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='desktop',
          debug=False,
          strip=None,
          upx=True,
          console=False )
coll = COLLECT(exe,
a.binaries + [('libwx_osx_cocoau-3.0.0.0.0.dylib',
               '/Users/jkitzes/anaconda/pkgs/wxpython-3.0-py27_0/lib/libwx_osx_cocoau-3.0.0.0.0.dylib',
               'BINARY')],
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='desktop')
app = BUNDLE(coll,
             name='desktop.app',
             icon=None)

