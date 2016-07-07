begin

  ; Create the file
  fname = "cvmix.nc"
  system("rm -f "+fname)
  fout = addfile(fname, "c")

  setfileoption(fout, "DefineMode", True)
  fAtt = True
  fAtt@title = "Standardized Output for CVMix Tests"
  fAtt@Conventions = "CMOR / CMIP6"
  fAtt@creation_date = systemfunc("date")
  fileattdef(fout, fAtt)

  ; Dimensions
  ; ---------
  ; time  = number of time steps (unlimited)
  ; nlev  = number of levels
  ; nface = number of interfaces = nlev + 1
  dimNames = (/"time", "nlev", "nface"/)
  dimSizes = (/    -1,    400,     401/)
  dimUnlim = (/  True,  False,   False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  ; Variables
  ; ---------
  ; time
  filevardef(fout,  "time", "double",  "time")
  t = new((/721/), "double")
  t = fspan(0, 864000, 721)
  t@units     = "seconds"
  t@long_name = "Elapsed Time"
  filevarattdef(fout, "time", t)

  filevardef(fout,   "zt", "double",   (/"time", "nlev"/))
  zt = new((/721,400/), "double")
  do n=0,720
    zt(n,:) = fspan(0.5d, 399.5d, 400)
  end do
  zt@units     = "m"
  zt@long_name = "Cell Center Depth"
  filevarattdef(fout, "zt", zt)

  filevardef(fout, "zw", "double", (/"time", "nface"/))
  zw = new((/721,401/), "double")
  do n=0,720
    zw(n,:) = fspan(0, 400, 401)
  end do
  zw@units     = "m"
  zw@long_name = "Cell Interface Depth"
  filevarattdef(fout, "zw", zw)

  new1d       = new((/721/),     "double")
  new2d_cntr  = new((/721,400/), "double")
  new2d_iface = new((/721,401/), "double")

  filevardef(fout, "blot", "double", "time")
  new1d@units     = "m"
  new1d@long_name = "Boundary Layer Ocean Thickness"
  filevarattdef(fout, "blot", new1d)

  ; T & S
  filevardef(fout, "thetao", "double", (/"time", "nlev"/))
  new2d_cntr@units     = "deg C"
  new2d_cntr@long_name = "Sea Water Potential Temperature"
  filevarattdef(fout, "thetao", new2d_cntr)

  filevardef(fout, "so", "double", (/"time", "nlev"/))
  new2d_cntr@units     = "Absolute Salinity (g/kg)"
  new2d_cntr@long_name = "Sea Water Salinity"
  filevarattdef(fout, "so", new2d_cntr)

  ; VELOCITY
  filevardef(fout, "uo", "double", (/"time", "nlev"/))
  new2d_cntr@units     = "m/s"
  new2d_cntr@long_name = "Sea Water x-Velocity"
  filevarattdef(fout, "uo", new2d_cntr)

  filevardef(fout, "vo", "double", (/"time", "nlev"/))
  new2d_cntr@units     = "m/s"
  new2d_cntr@long_name = "Sea Water y-Velocity"
  filevarattdef(fout, "vo", new2d_cntr)

  ; DIFFUSIVITY
  filevardef(fout, "difvho", "double", (/"time", "nface"/))
  new2d_cntr@units     = "m^2/s"
  new2d_cntr@long_name = "Ocean Vertical Heat Diffusivity"
  filevarattdef(fout, "difvho", new2d_cntr)

  filevardef(fout, "difvso", "double", (/"time", "nface"/))
  new2d_cntr@units     = "m^2/s"
  new2d_cntr@long_name = "Ocean Vertical Salt Diffusivity"
  filevarattdef(fout, "difvso", new2d_cntr)

  filevardef(fout, "difvmo", "double", (/"time", "nface"/))
  new2d_cntr@units     = "m^2/s"
  new2d_cntr@long_name = "Ocean Vertical Momentum Diffusivity"
  filevarattdef(fout, "difvmo", new2d_cntr)

  ; NON-LOCAL TERMS
  filevardef(fout, "difvhonl", "double", (/"time", "nface"/))
  new2d_cntr@units     = "unitless"
  new2d_cntr@long_name = "Non-Local Component of Ocean Vertical Heat Diffusivity"
  filevarattdef(fout, "difvhonl", new2d_cntr)

  filevardef(fout, "difvsonl", "double", (/"time", "nface"/))
  new2d_cntr@units     = "unitless"
  new2d_cntr@long_name = "Non-Local Component of Ocean Vertical Salt Diffusivity"
  filevarattdef(fout, "difvsonl", new2d_cntr)

  filevardef(fout, "difvmonl", "double", (/"time", "nface"/))
  new2d_cntr@units     = "unitless"
  new2d_cntr@long_name = "Non-Local Component of Ocean Vertical Momentum Diffusivity"
  filevarattdef(fout, "difvmonl", new2d_cntr)

  setfileoption(fout, "DefineMode", False)
  fout->time = (/t/)
  fout->zt   = (/zt/)
  fout->zw   = (/zw/)


end
