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
  shapefile = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/ganymed_model_new_shape_h.plt';
  trefdir = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/KRC/Ganymed/'
  outdir = '/Users/zouxd/Documents/GitHub/Thermal_fit/Ganymed/temp_maps2/'
  subelat = 36.3
  subslat = 18.9
  delta_lon = 27.1
  shape_a = 0.4250
  shape_b = 0.4250
  shape_c = 0.310
  pxlscl = 0.1 ; pixel scale in mas

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


; --- ADD THIS FIX ---
; If the max index in 'tri' equals the number of vertices, it's 1-based.
if max(tri) eq (size(vert))[2] then begin
    print, 'Correcting 1-based indexing in triangle array...'
    tri = tri - 1
endif
; --------------------

;Tfiles = file_search(trefdir + '/*krc.fits')  ; Ceres temperature model for test
Tfiles = file_search(trefdir + 'tref_gamma*_emis*.fits')  ; Juno temperature models

; fits keys to propagrate to output
keys = ['ti', 'emiss', 'rho', 'c', 'p_orb', 'p_rot', 'a_skin', 'd_skin']


; -----------------------------------------------------------------
; MODIFICATION 1:
; Longitude Logic (Manual Sequence)
; -----------------------------------------------------------------

subelon = [300.0, 332.0, 5.0, 56.0, 88.0, 127.0]
nlon = n_elements(subelon)

subelat = replicate(subelat, nlon)
subslon = (subelon + delta_lon) mod 360
subslat = replicate(subslat, nlon)
rh = replicate(rh_val, nlon)
delta = replicate(delta_val, nlon)

print, 'Output Longitudes will be: ', subelon
;start_lon = 27.0
; Calculate how many steps needed to cover 360 degrees roughly
;nlon = ceil(360.0 / dlon)

; Create the longitude array
;subelon = start_lon + findgen(nlon) * dlon
; Ensure it wraps around 360 (optional, keeps values 0-360)
;subelon = subelon mod 360.0

; Create arrays for geometry matching the size of subelon
;subelat = replicate(subelat, nlon)
;subslon = (subelon + delta_lon) mod 360
;subslat = replicate(subslat, nlon)
;rh = replicate(rh_val, nlon)
;delta = replicate(delta_val, nlon)

;print, 'Output Longitudes will be: ', subelon


for j=0, n_elements(Tfiles) - 1 do begin

  fname = (strsplit(Tfiles[j], '/', /ext))[-1]
  print, '--------------------------------------------------'
  print, 'Processing file: ', fname

  ; 1. Load reference temperature array (Main Data)
  tref = readfits(Tfiles[j], hdr, /silent)

  ; 2. Load LST (Time) & Reform
  lst = float(readfits(Tfiles[j], ext=1, /silent))
  lst = reform(lst, n_elements(lst))

  ; --- SANITY CHECK: Skip Corrupted Files ---
  ; If LST is huge (e.g., 8316 instead of ~48), the file is broken.
  if (n_elements(lst) gt 1000) then begin
    print, '>>> ERROR: Skipping ' + fname
    print, '    REASON: Corrupted Time Extension (Size: ' + strtrim(n_elements(lst),2) + ').'
    continue ; Jump to the next file immediately
  endif
  ; ------------------------------------------

  ; 3. Load ZZ (Depth) & Lats - Only proceed if file is good
  zz = float(readfits(Tfiles[j], ext=2, /silent))
  zz = reform(zz, n_elements(zz))

  lats = float(readfits(Tfiles[j], ext=3, /silent))
  lats = reform(lats, n_elements(lats))

  ; Calculate hour_angle
  hour_angle = ((lst - 12.0) * 15.0 + 360.0) mod 360.0
  time_arg = [hour_angle, 360.0]

  ; Prepare output directory
  dirname = outdir + strjoin((strsplit(file_basename(Tfiles[j]), '.', /ext))[0:-2], '.')
  if not file_test(dirname, /dir) then file_mkdir, dirname

  ; Loop throughout longitudes
  for i=0, n_elements(subslon)-1 do begin

    print, '  Sub-Earth longitude ', subelon[i]

    ; Calculate observing geometry
    sunpos = vect_match_view(rd2xyz([subslon[i], subslat[i]]), subelat[i], subelon[i], 0) * rh[i] * 1.496e8
    vert1 = vect_match_view(vert, subelat[i], subelon[i], 0)
    res = pxlscl / 206265000 * delta[i] * 1.496e8

    mesh_geomap, vert1, tri, sunpos, delta[i] * 1.496e8, imap, emap, amap, mask, pltmap, xres=res, yres=res, xs=ys, ys=ys

    ; Save plate map file
    suffix = '_' + string(round(subelon[i]), format='(i3.3)') + '.fits'
    outfile = dirname + '/platemap' + suffix

    mkhdr, hdr1, pltmap, /ext
    fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
    writefits, outfile, pltmap, hdr1
    writefits, outfile, emap, /append
    writefits, outfile, imap, /append
    writefits, outfile, amap, /append

    ; Calculate tpm
    temp_list = tpm_sph2plt(vert, tri, tref, subslon[i], time=time_arg, lats=lats)

    ; Temperature image mapping
    tempmap = tpm_mapping(temp_list, pltmap, depth=zz, zz=zz)

    ; Save simulation images
    outfile = dirname + '/tempmap' + suffix
    mkhdr, hdr1, tempmap, /ext
    fxaddpar, hdr1, 'bunit', 'K'
    for k=0, n_elements(keys) - 1 do begin
      val = sxpar(hdr, keys[k], comment=com)
      fxaddpar, hdr1, keys[k], val, com
    endfor
    fxaddpar, hdr1, 'long', subelon[i], 'longitude (deg)'
    writefits, outfile, tempmap, hdr1

    ; Append extensions
    mkhdr, hdr1, zz, /ima
    fxaddpar, hdr1, 'bunit', 'm'
    writefits, outfile, zz, hdr1, /append

    mkhdr, hdr1, emap, /ima
    fxaddpar, hdr1, 'bunit', 'deg'
    writefits, outfile, emap, hdr1, /append

  endfor

endfor

end
