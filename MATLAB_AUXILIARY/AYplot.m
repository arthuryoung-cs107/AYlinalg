classdef AYplot
    properties (Constant)
        fullpg_w = 1084;
        height1 = 200;
        margins_w = floor(AYplot.fullpg_w*0.8);

        posdim_full = [1 1 1728 1000];
        posdim_toprow = [1 551 1728 460];
        posdim_bottomrow = [1 1 1728 460];
        posdim_topleft=[1 551 576 460];
        posdim_topmiddle=[577 551 576 460];
        posdim_topright=[1153 551 576 460];
        posdim_movie = [10 30 500 500];

        pos0 = [1, 1];
        pos1 = [1, 750];
        pos2 = [1, 420];
        pos3 = [862, 750];
        pos4 = [862, 420];
        pos5 = [1, 442];
        pos6 = [862, 442];
        pos7 = [862, 1];

        size1 = [AYplot.margins_w, AYplot.height1];
        size2 = [AYplot.margins_w, 2*AYplot.height1];
        size3 = [AYplot.margins_w, 320];
        size4 = [AYplot.margins_w, 0.9*2*AYplot.height1];
        size5 = [AYplot.margins_w, 423];
        size6 = [400, 200];

        view_mat = [45, 45; 1, 0; 0, 90; 90, 0 ; 45, 0; 70, 10; -20, 10; -220, 10];

        plot_specs_def = struct('specs','o','color',[0.1,0.444,0.244],'lw',1,'ms',5);

        save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
        save_type = 'png';
    end
    methods
        function obj = AYplot()

        end
    end
    methods (Static)
        function AYfig_out = make_AYfig(name_,posdim_,make_ax_)
            if (nargin<3)
                make_ax=false;
            else
                make_ax=make_ax_;
            end
            AYfig_out = AYfig(AYfig.specs_gen(name_,posdim_),make_ax);
        end
        function plot3_scaled_vectors(axs_,p_,v_,sp_,normalization_flag_)
            if (nargin==4)
                normalization_flag=false;
            else
                normalization_flag=normalization_flag_;
            end

            [xl yl zl] = deal(xlim(axs_),ylim(axs_),zlim(axs_));
            lim_diff = [diff(xl) diff(yl) diff(zl)];
            max_diff = max(lim_diff);
            s = lim_diff/max_diff;
            vs = v_.*s;
            if (normalization_flag)
                p = p_+((vs./(sqrt(sum(vs.*vs,2)))).*s);
            else
                p = p_+(vs.*s);
            end

            plot3(axs_, [p_(:, 1), p(:, 1)]', [p_(:, 2), p(:, 2)]', [p_(:, 3), p(:, 3)]',sp_.specs, 'Color',sp_.color,'LineWidth',sp_.lw,'MarkerSize',sp_.ms);
        end
    end
end
