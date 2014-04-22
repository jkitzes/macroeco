# -*- mode: python -*-
a = Analysis(['mecodesktop.py'],
             pathex=['/Users/jkitzes/Projects/macroeco'],
             hiddenimports=['scipy.special._ufuncs_cxx'],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='mecodesktop',
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
               name='mecodesktop')
app = BUNDLE(coll,
             name='Macroeco Desktop.app',
             icon='icon.icns')

