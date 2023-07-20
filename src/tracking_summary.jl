struct TrackingSummary{T <: AbstractGNSS}
    constellation :: T
    prn::Int
    channel_id::Int
    acq_delay_samples::Float64
    acq_doppler_hz::Float64
    acq_samplestamp_samples::UInt64
    acq_doppler_step::UInt32
    flag_valid_acquisition::Bool
    fs::Float64
    prompt::ComplexF64
    cn0_db_hz::Float64
    carrier_doppler_hz::Float64
    carrier_phase_rads::Float64
    code_phase_samples::Float64
    tracking_sample_counter::UInt64
    flag_valid_symbol_output::Bool
    correlation_length_ms::Int32
end

function TrackingSummary(trk, channel_id, samplerate)
    return TrackingSummary(
        get_state(trk).system,
        get_state(trk).prn,
        channel_id,
        0.0,
        0.0,
        UInt64(0),
        UInt32(0),
        true,
        samplerate,
        get_prompt(trk),
        ustrip(trk.cn0),
        ustrip(get_carrier_doppler(trk)),
        get_carrier_phase(trk),
        get_code_phase(trk),
        UInt64(0),
        true,
        Int32(1)
    )
end

function Base.show(io::IO, ::MIME"text/plain", trk_channels::Dict{Int64, GNSSNavigation.TrackingSummary})
    header = ["Channel"; "PRN"; "CN0"; "Carrier doppler (Hz)"; "Code phase (samples)"]
    data = Matrix{Any}(undef, length(trk_channels),length(header))

    for (idx,prn_val) in enumerate(trk_channels)
        data[idx,1] = idx
        data[idx,2] = prn_val[1]
        data[idx,3] = prn_val[2].cn0_db_hz
        data[idx,4] = prn_val[2].carrier_doppler_hz
        data[idx,5] = prn_val[2].code_phase_samples
    end
    hl_good = Highlighter((data,i,j)->(j==3) &&(data[i,j] > 42),crayon"green")
    hl_bad = Highlighter((data,i,j)->(j==3) &&(data[i,j] < 42),crayon"red")
    
    pretty_table(io,data,header=header,highlighters=(hl_good,hl_bad))
end
