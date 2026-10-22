using Test
include(joinpath(@__DIR__, "..", "example", "cycling_csv_time.jl"))

@testset "CSV mapping uses cycle and phase identity, never time resets" begin
    phase(start, duration) = (; t_start=start, duration,
        solve_result=Dict("time [s]"=>[0.0, duration]))
    # Rest1 has no snapshot; cycle 2 has no rest1 at all. Equal local stamps
    # must still map to their real phase and cumulative cycle origin.
    c1=(cycle_idx=1,discharge=phase(0.0,1.0),rest1=phase(1.0,1.0),
        charge=phase(2.0,1.0),rest2=phase(3.0,1.0))
    c2=(cycle_idx=2,discharge=phase(4.0,1.0),rest1=nothing,
        charge=phase(5.0,1.0),rest2=phase(6.0,1.0))
    snaps=[(cycle=2,phase="rest2",time_s=1.0),
        (cycle=1,phase="rest2",time_s=1.0),
        (cycle=1,phase="discharge",time_s=1.0),
        (cycle=2,phase="charge",time_s=1.0)]
    r=(cycle_results=[c1,c2],czm_snapshots=snaps)
    mapped=cycling_csv_snapshots(r)
    @test [x.time_s for x in mapped]==[1.0,4.0,6.0,7.0]
    @test [x.phase for x in mapped]==["discharge","rest2","charge","rest2"]
    @test mapped[1].snapshot===snaps[3]
    @test cycling_csv_times(c2.charge,2)==[5.0,6.0]
    @test_throws DimensionMismatch cycling_csv_times(c2.charge,3)
    @test_throws KeyError cycling_csv_times((t_start=0.0,solve_result=Dict()),2)
    @test_throws KeyError cycling_csv_snapshots((;r...,czm_snapshots=[(cycle=2,phase="rest1",time_s=1.0)]))
    @test_throws KeyError cycling_csv_snapshots((;r...,czm_snapshots=[(cycle=1,phase="PHASE_REST",time_s=1.0)]))
    @test_throws ErrorException cycling_csv_snapshots((;r...,czm_snapshots=[(cycle=1,phase="rest2",time_s=2.0)]))
    @test isempty(cycling_csv_snapshots((;r...,czm_snapshots=snaps[1:0])))
end
