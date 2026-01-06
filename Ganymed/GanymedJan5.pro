; Script to simulate temperature image cubes of Didymos/Dimorphos based on the
; pre-calculated T_ref from T.N. Titus
;
; Models are pre-processed to break down to individual files each containing a
; temperature model in the format fully complying with the requirement of
; 'tpm_sph2plt.pro'.
;
; Input temperatur emodel are saved in subdirectory 'KRC/'.  Each individual
; model file has a name like 'tref_gamma{:d}_emis{:.1f}.fits'.format(ti, emis).


; switch for Juno shape model or triaxial ellipsoid
shape = 'gany'   ; use Juno shape model

if shape eq 'gany' then begin
  shapefile = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/Ganymed_261shape_tmp.plt';
  trefdir = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/KRC/Ganymed/'
  outdir = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/temp_maps2/'
  subelat = 36.3
  subslat = 18.9
  delta_lon = 27.1
  shape_a = 0.4250
  shape_b = 0.4250
  shape_c = 0.310
  pxlscl = 2.8 ; pixel scale in mas

endif

; longitude step size for full lightcurve calculation
dlon = 33.0      ; longitudinal step size for rotation (Use float)

; latitudes
rh_val = 1.31 ; solar distance (renamed to avoid confusion with array)
delta_val = 0.39 ; range from Earth (renamed to avoid confusion with array)

; simulated image parameters

xs = 256   ; x size of image
ys = 256   ; y size of image


;-------------------------------------------

; load shape model
readshape_triplate, shapefile, vert, tri
print, 'Dimensions of vert, tri: ', size(vert), size(tri)

;Tfiles = file_search(trefdir + '/*krc.fits')  ; Ceres temperature model for test
Tfiles = file_search(trefdir + 'tref_gamma*_emis*.fits')  ; Juno temperature models

; fits keys to propagrate to output
keys = ['ti', 'emiss', 'rho', 'c', 'p_orb', 'p_rot', 'a_skin', 'd_skin']


; -----------------------------------------------------------------
; MODIFICATION 1: Longitude Logic
; Start at 27 deg, Step 33 deg.
; -----------------------------------------------------------------
start_lon = 27.0
; Calculate how many steps needed to cover 360 degrees roughly
nlon = ceil(360.0 / dlon) 

; Create the longitude array
subelon = start_lon + findgen(nlon) * dlon
; Ensure it wraps around 360 (optional, keeps values 0-360)
subelon = subelon mod 360.0

; Create arrays for geometry matching the size of subelon
subelat = replicate(subelat, nlon)
subslon = (subelon + delta_lon) mod 360
subslat = replicate(subslat, nlon)
rh = replicate(rh_val, nlon)
delta = replicate(delta_val, nlon)

print, 'Output Longitudes will be: ', subelon

for j=0, n_elements(Tfiles) - 1 do begin

  print, 'processing file: '+(strsplit(Tfiles[j], '/', /ext))[-1]

  ; load reference temperature array
  tref = readfits(Tfiles[j], hdr, /silent)
  
  ; -----------------------------------------------------------------
  ; MODIFICATION 2: Fix Dimensions (Debug)
  ; Use REFORM() to ensure arrays are 1D vectors, not (1, N) or (N, 1)
  ; -----------------------------------------------------------------
  lst = readfits(Tfiles[j], ext=1, /silent)
  lst = reform(lst) ; <--- CRITICAL FIX for concatenation error
  
  hour_angle = ((lst - 12) * 15 + 360) mod 360  ; hour angle from local noon in degrees
  
  zz = readfits(Tfiles[j], ext=2, /silent)
  zz = reform(zz)   ; <--- Good practice to reform this too
  
  lats = readfits(Tfiles[j], ext=3, /silent)
  lats = reform(lats) ; <--- Good practice to reform this too

  ; prepare output directory
  dirname = outdir + strjoin((strsplit(file_basename(Tfiles[j]), '.', /ext))[0:-2], '.')
  if not file_test(dirname, /dir) then file_mkdir, dirname

  ; loop throughout longitudes
  for i=0, n_elements(subslon)-1 do begin 

    print, 'Sub-Earth longitude ', subelon[i]
    
    ; calculate observing geometry
    sunpos = vect_match_view(rd2xyz([subslon[i], subslat[i]]), subelat[i], subelon[i], 0) * rh[i] * 1.496e8
    vert1 = vect_match_view(vert, subelat[i], subelon[i], 0)
    res = pxlscl / 206265000 * delta[i] * 1.496e8   ; image resolution in km/s
    
    mesh_geomap, vert1, tri, sunpos, delta[i] * 1.496e8, imap, emap, amap, mask, pltmap, xres=res, yres=res, xs=ys, ys=ys
    
    ; save plate map file
    suffix = '_' + string(round(subelon[i]), format='(i3.3)') + '.fits'

    ; save plate files
    outfile = dirname + '/platemap' + suffix
    mkhdr, hdr1, pltmap, /ext
    fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
    writefits, outfile, pltmap, hdr1
    writefits, outfile, emap, /append
    writefits, outfile, imap, /append
    writefits, outfile, amap, /append

    ; calculate tpm for plate shape model
    ; Because hour_angle is now reformed (1D), [hour_angle, 360] will work correctly
    temp_list = tpm_sph2plt(vert, tri, tref, subslon[i], time=[hour_angle, 360], lats=lats)
    
    ; temperature image
    tempmap = tpm_mapping(temp_list, pltmap, depth=zz, zz=zz)

    ; save simulation images
    outfile = dirname + '/tempmap' + suffix
    ; prepare primary extension
    mkhdr, hdr1, tempmap, /ext
    fxaddpar, hdr1, 'bunit', 'K'
    for k=0, n_elements(keys) - 1 do begin
      val = sxpar(hdr, keys[k], comment=com)
      fxaddpar, hdr1, keys[k], val, com
    endfor
    fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
    writefits, outfile, tempmap, hdr1
    ; first extension: depth
    mkhdr, hdr1, zz, /ima
    fxaddpar, hdr1, 'bunit', 'm'
    writefits, outfile, zz, hdr1, /append
    ; second extension: emission angle
    mkhdr, hdr1, emap, /ima
    fxaddpar, hdr1, 'bunit', 'deg'
    writefits, outfile, emap, hdr1, /append

  endfor

endfor

end
