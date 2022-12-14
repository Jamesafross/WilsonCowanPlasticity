include("functions/WC_InitSetup.jl")

global WCS = setup()

WCS.modelparams.ρ=0.3
WCS.solver_opts.tWindow = 400
#WCS.solver_opts.use_noise_process = false
#WCS.solver_opts.N = NoiseWrapper(sol.W)
@time wilsoncowan_windows_run()

modelFC_nostim = getModelFC(WCS.bold_out,1,1) 

save_sol = zeros(size(sol[1:5*140,:]))
save_sol .= sol[1:5*140,:]

global WCS = setup()


WCS.modelparams.ρ=0.3
WCS.solver_opts.tWindow = 4000
WCS.stimparams.stim = true
WCS.stimparams.Tstim = [500,530]
WCS.stimparams.stimStr = -5
WCS.solver_opts.N = NoiseWrapper(sol.W)
WCS.solver_opts.run = 2


@time wilsoncowan_windows_run()

modelFC_stim = getModelFC(WCS.bold_out,1,1) 

diff = abs.(modelFC_nostim.^2 .- modelFC_stim.^2)
diffpercent = (diff./((modelFC_nostim.^2)))*100


function plotFM(window)
    p1=heatmap(modelFC_nostim[:,:,window].^2)
    p2=heatmap(modelFC_stim[:,:,window].^2)
    p3 = heatmap(abs.(modelFC_nostim[:,:,window].^2 .- modelFC_stim[:,:,window].^2))
    return plot(p1,p2,p3)
end
