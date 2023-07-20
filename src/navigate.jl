#= function navigate(fup, n, acq_results; intermediate_freq=0, dsp_callback=nothing, userdata=nothing)
    #resampler_ctx_ch0 = FIRFilter(7//4)
    fp = read(fup, 10u"ms")
    #fp = Kea.resample(fp[:,1], resampler_ctx_ch0)
    trk_results = similar(acq_results, Tracking.TrackingResults)
    trk_states = similar(acq_results, Tracking.TrackingState)
    decoder_states = similar(acq_results, typeof(GPSL1DecoderState(1)))
    subctr = 0
    subctr1 = 0 
    Hz = u"Hz"
    for i in eachindex(acq_results)
      prn = acq_results[i].prn
      decoder_states[i] = GPSL1DecoderState(prn)
      trk_states[i] = TrackingState(prn, GPSL1(), acq_results[i].carrier_doppler, acq_results[i].code_phase)
    end
    #init states
    Threads.@threads for i in eachindex(trk_results)
      trk_results[i] = track(ComplexF64.(collect(fp[:,1])), trk_states[i], fp.samplerate*Hz;intermediate_frequency=intermediate_freq*Hz)
    end
    filtered_prompt = Vector{ComplexF64}()
    clock_drifts = Vector{Float64}()
    trk_result_out = Vector{Vector{Tracking.TrackingResults}}()
    pvt_sol = Vector{PVTSolution}()
    file_idx = Vector{Int}()
    ret = Vector{Any}()

    for kk in 1:n
      fp = read(fup, 10u"ms")
      #fp = Kea.resample(fp[:,1], resampler_ctx_ch0)
      #println("Iteration $kk")
      #samples = ComplexF64.(fp[:,1].data);
      rate = fp.samplerate*Hz;
      #dat = Ref(fp[:,1].data)
      Threads.@threads for i in eachindex(trk_results)
        dat = fp[:,1].data
        trk_results[i] = track(dat, get_state(trk_results[i]), rate;intermediate_frequency=intermediate_freq*Hz);
        #append!(filtered_prompt,trk_res.filtered_prompt)
        #println(trk_results[i].cn0)
        decoder_states[i] = decode(decoder_states[i], get_bits(trk_results[i]), get_num_bits(trk_results[i]))
      end

      if subctr1 % 50 == 0
          for i in 1:6
            #@info "PRN $(trk_results[i].state.prn) $(trk_results[i].cn0)"
          end
      end
      subctr1 = subctr1 + 1

      qq = SatelliteState.(decoder_states, trk_results);
      #println(qq)
      nav = calc_pvt(qq);
      if !(isnothing(nav.time))
        lla = get_LLA(nav)
        if subctr == 0
          @info "First nav soln at $(nav.time) is lat: $(lla.lat) lon: $(lla.lon) alt: $(lla.alt) clk drift: $(nav.relative_clock_drift)"
        else
          if subctr % 50 == 0
            @info "Position at time $(nav.time) is lat: $(lla.lat) lon: $(lla.lon) alt: $(lla.alt) clk drift: $(nav.relative_clock_drift)"
            @info "Vx $(nav.velocity.x) m/s Vy $(nav.velocity.y) m/s Vz $(nav.velocity.z) m/s "
          end
        end
       # append!(trk_result_out, [trk_results])
        #append!(pvt_sol, [nav])
        #append!(file_idx,[kk])
        #append!(clock_drifts,nav.relative_clock_drift)

        a = []
        for i in 1:length(trk_results)
            append!(a,[TrackingSummary(trk_results[i],i,fp.samplerate)])
        end

        append!(ret, ((kk, a, nav),))
        subctr = subctr +1
      end
    end
  #return pvt_sol,trk_result_out, file_idx
  return ret
end =#

function navigate(dat, n, acq_results, samplerate; intermediate_freq=0, dsp_callback=nothing, userdata=nothing, silent=false) where T <: Number
  #resampler_ctx_ch0 = FIRFilter(7//4)
  #fp = Kea.resample(fp[:,1], resampler_ctx_ch0)
  trk_results = similar(acq_results, Tracking.TrackingResults)
  trk_states = similar(acq_results, Tracking.TrackingState)
  decoder_states = similar(acq_results, typeof(GPSL1DecoderState(1)))
  subctr = 0
  subctr1 = 0 
  Hz = u"Hz"
  ms_10 = Int(fld(samplerate,100))
  for i in eachindex(acq_results)
    prn = acq_results[i].prn
    decoder_states[i] = GPSL1DecoderState(prn)
    trk_states[i] = TrackingState(prn, GPSL1(), acq_results[i].carrier_doppler, acq_results[i].code_phase)
  end
  #init states
  Threads.@threads for i in eachindex(trk_results)
    @views trk_results[i] = track(ComplexF32.(dat[1:ms_10]), trk_states[i], samplerate*Hz;intermediate_frequency=intermediate_freq*Hz)
  end
  filtered_prompt = Vector{ComplexF64}()
  clock_drifts = Vector{Float64}()
  trk_result_out = Vector{Vector{Tracking.TrackingResults}}()
  pvt_sol = Vector{PVTSolution}()
  file_idx = Vector{Int}()
  decodes = []
  ret = Vector{Any}()
  header = ["Channel"; "PRN"; "CN0"; "Carrier doppler (Hz)"; "Code phase (samples)"; "Nav msg"]
  data = Matrix{Any}(undef, length(trk_results),length(header))


  for kk in 2:n
    #fp = Kea.resample(fp[:,1], resampler_ctx_ch0)
    #println("Iteration $kk")
    #samples = ComplexF64.(fp[:,1].data);
    rate = samplerate*Hz;
    #dat = Ref(fp[:,1].data)
    Threads.@threads for i in eachindex(trk_results)
      @views trk_results[i] = track(dat[((ms_10*kk)+1):(ms_10*(kk+1))], get_state(trk_results[i]), rate;intermediate_frequency=intermediate_freq*Hz);
      #append!(filtered_prompt,trk_res.filtered_prompt)
      #println(trk_results[i].cn0)
      decoder_states[i] = decode(decoder_states[i], get_bits(trk_results[i]), get_num_bits(trk_results[i]))
    end
    

    if subctr1 % 50 == 0
        for trk in trk_results
          #if !silent println("PRN $(trk.state.prn) $(trk.cn0)") end
        end
    end
    subctr1 = subctr1 + 1

    qq = SatelliteState.(decoder_states, trk_results);
    #println(qq)
    nav = calc_pvt(qq);

    a = Dict{Int,TrackingSummary}()
    for i in 1:length(trk_results)
        #append!(a,[TrackingSummary(trk_results[i],i,Float64(samplerate))])
        a[get_state(trk_results[i]).prn] = TrackingSummary(trk_results[i],i,Float64(samplerate))
    end

    if !(isnothing(nav.time))
      lla = get_LLA(nav)
      if subctr == 0
        #if !silent println("First nav soln at $(nav.time) is lat: $(lla.lat) lon: $(lla.lon) alt: $(lla.alt) clk drift: $(nav.relative_clock_drift)") end
      else
        if subctr % 50 == 0
 #=          if !silent
            println("Position at time $(nav.time) is lat: $(lla.lat) lon: $(lla.lon) alt: $(lla.alt) clk drift: $(nav.relative_clock_drift)")
            println("Vx $(nav.velocity.x) m/s Vy $(nav.velocity.y) m/s Vz $(nav.velocity.z) m/s ")
          end =#
        end
      end
     # append!(trk_result_out, [trk_results])
      #append!(pvt_sol, [nav])
      #append!(file_idx,[(ms_10*kk)+1])
      #append!(clock_drifts,nav.relative_clock_drift)

 #=      a = Dict{Int,TrackingSummary}()
      for i in 1:length(trk_results)
          #append!(a,[TrackingSummary(trk_results[i],i,Float64(samplerate))])
          a[get_state(trk_results[i]).prn] = TrackingSummary(trk_results[i],i,Float64(samplerate))
      end
 =#
      #append!(ret, (((ms_10*kk)+1, a, nav, deepcopy(decoder_states)),))
      subctr = subctr +1
    end

    append!(ret, (((ms_10*kk)+1, a, nav, ()),))

    if !silent
      for (idx,prn_val) in enumerate(a)
        data[idx,1] = idx
        data[idx,2] = prn_val[1]
        data[idx,3] = prn_val[2].cn0_db_hz
        data[idx,4] = prn_val[2].carrier_doppler_hz
        data[idx,5] = prn_val[2].code_phase_samples
        #n^2 search... todo fix this
        for decoder_state in decoder_states
          if decoder_state.prn == prn_val[1]
            navmsg = [" ","  ", " ","  ", " "]
            if GNSSDecoder.is_subframe1_decoded(decoder_state.raw_data)
              navmsg[1] = "1"
            end
            if GNSSDecoder.is_subframe2_decoded(decoder_state.raw_data)
              navmsg[3] = "2"
            end
            if GNSSDecoder.is_subframe3_decoded(decoder_state.raw_data)
              navmsg[5] = "3"
            end
            data[idx,6] = prod(navmsg)
          end
        end
  
      end
      hl_good = Highlighter((data,i,j)->(j==3) &&(data[i,j] > 42),crayon"green")
      hl_bad = Highlighter((data,i,j)->(j==3) &&(data[i,j] < 42),crayon"red")

      nav_bad = Highlighter((data,i,j)->(j==6) &&(data[i,j] != "1  2  3"),crayon"red")
      nav_good = Highlighter((data,i,j)->(j==6) &&(data[i,j] == "1  2  3"),crayon"green")
      
      pretty_table(data,header=header,highlighters=(hl_good,hl_bad,nav_bad,nav_good),overwrite=false)
      if !(isnothing(nav.time))
        lla = get_LLA(nav)
        println("Position at time $(nav.time) is lat: $(lla.lat) lon: $(lla.lon) alt: $(lla.alt) clk drift: $(nav.relative_clock_drift)")
        println("Vx $(nav.velocity.x) m/s Vy $(nav.velocity.y) m/s Vz $(nav.velocity.z) m/s ")
      end
    end
  end
#return pvt_sol,trk_result_out, file_idx
return ret
end