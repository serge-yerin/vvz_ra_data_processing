Project structure                                                     
                                                                        
  idl_to_python/                                                        
  ├── indsearch/                                            
  │   ├── __init__.py                                                   
  │   ├── dispersion.py   — DM delay formula                            
  │   ├── reader.py       — binary file I/O                             
  │   ├── processor.py    — main dedispersion loop                      
  │   └── plot.py         — matplotlib DM-time visualization            
  ├── main.py             — CLI entry point                             
  └── requirements.txt                                                  
                                                                        
  Run it                                                                
                                                                        
  pip install -r requirements.txt
  python3 main.py observation.ucd 26.7                                  
  python3 main.py observation.ucd 26.7 --no-plot   # skip figure        

  Key translation decisions

  IDL construct: fltarr(wofsg, nofs) read via assoc
  Python equivalent: np.frombuffer(...).reshape((wofsg,nofs), order='F')
  Why: IDL stores arrays in Fortran (column-major) order; order='F'
  makes
    frame[i,j] match IDL's data[i,j]
  ────────────────────────────────────────
  IDL construct: total(arr, 1)
  Python equivalent: arr.sum(axis=0)
  Why: IDL dimensions are 1-indexed; dim 1 → axis 0
  ────────────────────────────────────────
  IDL construct: shift(arr, n)
  Python equivalent: np.roll(arr, n)
  Why: Identical circular-shift semantics and sign convention
  ────────────────────────────────────────
  IDL construct: writeu 2D array
  Python equivalent: out_slice.T.astype('<f4').tofile()
  Why: IDL writeu writes Fortran order (DM index varies fastest);
    transposing before tofile() reproduces that byte layout
  ────────────────────────────────────────
  IDL construct: arr[6:32] (inclusive)
  Python equivalent: arr[6:33]
  Why: IDL slice endpoints are inclusive; Python's are exclusive
  ────────────────────────────────────────
  IDL construct: ULONARR(128)
  Python equivalent: dtype='<u4'
  Why: Unsigned 32-bit little-endian

