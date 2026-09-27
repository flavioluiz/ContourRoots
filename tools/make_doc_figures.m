function make_doc_figures()
%MAKE_DOC_FIGURES Create the figures of docs/assets.
%   Draws the figures used only by the Markdown documentation and copies
%   selected study figures from output/ (run "buildtool examples" first).
%   The beam schematic is converted from TikZ when pdflatex and pdftocairo
%   are available.

    root = fileparts(fileparts(mfilename('fullpath')));
    assets = fullfile(root,'docs','assets');
    if ~isfolder(assets), mkdir(assets); end
    run(fullfile(root,'setup_contourroots.m'));
    old = get(groot,{'defaultFigureVisible','defaultAxesFontSize'});
    set(groot,'defaultFigureVisible','off','defaultAxesFontSize',12);
    restore = onCleanup(@() set(groot,{'defaultFigureVisible','defaultAxesFontSize'},old));

    overview(fullfile(assets,'readme_overview.png'));
    argument_principle(fullfile(assets,'argument_principle.png'));
    first_roots(fullfile(assets,'first_roots.png'));

    copies = {
        'time_delay/figures/critical_delays.png'
        'time_delay/figures/spectral_abscissa.png'
        'time_delay/figures/pole_maps.png'
        'time_delay/figures/pade_phase_universal.png'
        'time_delay/figures/pade_step_response.png'
        'time_delay/figures/p1_time_simulation.png'
        'time_delay/figures/p3_unstable_count.png'
        'time_delay/figures/pade_critical_phase.png'
        'distributed_systems/heat_and_characteristic.png'
        'distributed_systems/wave_duct_beam.png'
        'distributed_systems/beam_accumulation.png'
        'coupled_beam/pole_sweeps.png'
        'coupled_beam/trends_and_validation.png'
        'coupled_beam/coupled_modes.png'
        'time_response/time_response_validation.png'
        'time_response/time_response_aeroelasticity.png'
        'time_response/time_response_convergence.png'
        'continuous_wing/wing_root_locus.png'
        'continuous_wing/wing_convergence.png'
        'continuous_wing/wing_time_response.png'
        'hybrid_comparison/hybrid_comparison.png'
        'hybrid_comparison/hybrid_superposition.png'};
    for k = 1:numel(copies)
        src = fullfile(root,'output',copies{k});
        if ~isfile(src)
            error('make_doc_figures:MissingStudy', ...
                '%s not found. Run "buildtool examples" first.', src);
        end
        copyfile(src, assets);
    end
    aeroCopies = {
        'nasa_study.png', 'aeroelastic_nasa.png';
        'dlr_study.png', 'aeroelastic_dlr.png';
        'theodorsen_and_domain.png', 'aeroelastic_theodorsen.png';
        'nasa_parameter_study.png', 'aeroelastic_parameters.png'};
    for k = 1:size(aeroCopies,1)
        src = fullfile(root,'output','aeroelasticity',aeroCopies{k,1});
        if ~isfile(src)
            error('make_doc_figures:MissingStudy', ...
                '%s not found. Run "buildtool examples" first.', src);
        end
        copyfile(src,fullfile(assets,aeroCopies{k,2}));
    end
    schematic(root, assets);
    fprintf('Documentation figures written to %s\n', assets);
end

function overview(file)
    f = figure('Position',[100 100 1100 430],'Color','w');
    t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t); hold(ax,'on');
    F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
    region = [-8 2 -20 20];
    r = croots(F,region,'AssumeAnalytic',true);
    frame = region([1 2 2 1 1]) + 1i*region([3 3 4 4 3]);
    patch(ax,[0 2 2 0],[-20 -20 20 20],[1 .92 .9],'EdgeColor','none');
    plot(ax,real(frame),imag(frame),'--','Color',[.55 .55 .55]);
    xline(ax,0,':'); yline(ax,0,':');
    plot(ax,real(r),imag(r),'x','MarkerSize',11,'LineWidth',2,'Color',[0 .447 .741]);
    text(ax,0.15,17,'unstable','Color',[.7 .2 .1]);
    xlim(ax,[-8.5 2.5]); ylim(ax,[-21 21]); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Re(s)'); ylabel(ax,'Im(s)');
    title(ax,'croots: s^2+s+1+(2s+3)e^{-s} = 0');
    ax = nexttile(t);
    G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
    cpzmap(ax,G,[-1 1 -10 10],'AssumeAnalytic',true);
    title(ax,'cpzmap: sinh(s/2) / sinh(s)');
    exportgraphics(f,file,'Resolution',160); close(f);
end

function argument_principle(file)
    % Winding of F along the boundary of a rectangle equals the number of
    % zeros inside it.
    F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
    rect = [-2 1 -3 9];
    n = 4000; t = (0:n)/n;
    z = [rect(1)+diff(rect(1:2))*t+1i*rect(3), rect(2)+1i*(rect(3)+diff(rect(3:4))*t), ...
         rect(2)-diff(rect(1:2))*t+1i*rect(4), rect(1)+1i*(rect(4)-diff(rect(3:4))*t)];
    w = F(z);
    r = croots(F,[-8 2 -20 20],'AssumeAnalytic',true);
    inside = r(real(r)>rect(1) & real(r)<rect(2) & imag(r)>rect(3) & imag(r)<rect(4));
    f = figure('Position',[100 100 1250 400],'Color','w');
    t3 = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    ax = nexttile(t3); hold(ax,'on');
    plot(ax,real(z),imag(z),'-','Color',[.2 .2 .2],'LineWidth',1.5);
    quiver(ax,0,rect(3),0.5,0,0,'Color',[.2 .2 .2],'MaxHeadSize',2,'LineWidth',1.5);
    plot(ax,real(r),imag(r),'x','Color',[.65 .65 .65],'MarkerSize',9,'LineWidth',1.5);
    plot(ax,real(inside),imag(inside),'x','Color',[0 .447 .741],'MarkerSize',11,'LineWidth',2);
    xlim(ax,[-3 2]); ylim(ax,[-4 10]); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Re(s)'); ylabel(ax,'Im(s)');
    title(ax,sprintf('(a) Contour C with %d roots inside',numel(inside)));
    ax = nexttile(t3); hold(ax,'on');
    rho = log1p(abs(w)); wc = rho.*exp(1i*angle(w));
    plot(ax,real(wc),imag(wc),'-','Color',[.85 .33 .1],'LineWidth',1.3);
    plot(ax,0,0,'k+','MarkerSize',12,'LineWidth',1.5);
    axis(ax,'equal'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Re F'); ylabel(ax,'Im F');
    title(ax,'(b) Image F(C), radius log(1+|F|)');
    ax = nexttile(t3);
    arc = [0 cumsum(abs(diff(z)))];
    phase = unwrap(angle(w));
    plot(ax,arc,(phase-phase(1))/(2*pi),'-','Color',[0 .447 .741],'LineWidth',1.6);
    yline(ax,numel(inside),'--',sprintf('%d turns',numel(inside)));
    grid(ax,'on'); box(ax,'on'); xlim(ax,[0 arc(end)]);
    xlabel(ax,'arc length along C'); ylabel(ax,'change of arg F / 2\pi');
    title(ax,'(c) Winding number = number of roots');
    exportgraphics(f,file,'Resolution',160); close(f);
end

function first_roots(file)
    f = figure('Position',[100 100 560 430],'Color','w');
    ax = axes(f); hold(ax,'on');
    F = @(s) s + 1 + 2*exp(-1.5*s);
    [r,~] = croots(F,[-6 1 -30 30],'AssumeAnalytic',true);
    patch(ax,[0 1 1 0],[-30 -30 30 30],[1 .92 .9],'EdgeColor','none');
    plot(ax,real(r),imag(r),'x','MarkerSize',10,'LineWidth',2,'Color',[0 .447 .741]);
    xline(ax,0,':'); grid(ax,'on'); box(ax,'on');
    xlim(ax,[-6 1]); ylim(ax,[-30 30]);
    xlabel(ax,'Re(s)'); ylabel(ax,'Im(s)');
    title(ax,'Roots of s + 1 + 2e^{-1.5s}: two unstable');
    exportgraphics(f,file,'Resolution',160); close(f);
end

function schematic(root, assets)
    % TikZ figures of the manual, converted to SVG for the Markdown docs.
    figures = {'beam_schematic_standalone','beam_schematic.svg';
               'typical_section_standalone','aeroelastic_section.svg';
               'continuous_wing_standalone','continuous_wing.svg';
               'hybrid_ports_standalone','hybrid_ports.svg'};
    for k = 1:size(figures,1)
        tex = fullfile(root,'manual','figures',[figures{k,1} '.tex']);
        build = tempname; mkdir(build);
        cleaner = onCleanup(@() rmdir(build,'s'));
        cmd = sprintf(['cd "%s" && pdflatex -interaction=nonstopmode -output-directory="%s" ' ...
            '"%s" >/dev/null && pdftocairo -svg "%s" "%s"'], fileparts(tex), build, tex, ...
            fullfile(build,[figures{k,1} '.pdf']), fullfile(assets,figures{k,2}));
        [status,out] = system(['PATH=$PATH:/Library/TeX/texbin:/opt/homebrew/bin:/usr/local/bin; ' cmd]);
        if status ~= 0
            warning('make_doc_figures:Schematic','%s not converted:\n%s',figures{k,2},out);
        end
        clear cleaner
    end
end
