module GNSSNavigation

using Acquisition
using PositionVelocityTime
using GNSSSignals
using Tracking
using Unitful
using PrecompileTools
using GNSSDecoder
using PrettyTables
using TrackingSummaryDefinition


include("navigate.jl")
export navigate

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    rate = 2.048e6 

    data_ci8 = rand(Complex{Int8},Int(round(rate*12)))
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        samplerate = rate*u"Hz"
        a = coarse_fine_acquire(GPSL1(),data_ci8[1:2048], samplerate, 1:32; interm_freq=0*u"Hz", max_doppler=40e3*u"Hz");
        navigate(data_ci8, 1000, a[1:6], rate;silent=true);

    end
end


end
