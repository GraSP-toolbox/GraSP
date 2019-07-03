% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2019)
% 
% benjamin.girault@ens-lyon.fr
% benjamin.girault@usc.edu
% 
% This software is a computer program whose purpose is to provide a Matlab
% / Octave toolbox for handling and displaying graph signals.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

%% Graph Translation example
%
% We show here how the graph translation is constructed and how it compares
% to two different time shift like operators, namely the generalized
% translations and the graph shift.

%% Example graph
%
% We use a graph of $100$ vertices randomly drawn in a 2D square of size
% 10x10 using the uniform distribution.

g = grasp_plane_rnd(100);
g.show_graph_options.color_map = flipud(colormap('hot'));
g.show_graph_options.show_colorbar = true;
grasp_show_graph(gca, g, 'show_edges', false, 'show_colorbar', false);

%%
% Add edges between vertices at distance less than $3$. Weight the edges
% by a Gaussian kernel $a_{ij}=exp(-\frac{d(i,j)^2}{2\sigma^2})$ with
% $\sigma^2=1/1.5$. Note that grasp_adjacency_thresh thresholds the weights
% on the edges, and not the distances.

g.A = grasp_adjacency_gaussian(g, 1/1.5);
g.A = grasp_adjacency_thresh(g, exp(-3 * 3 ^ 2));
g.A_layout = 0;
grasp_show_graph(gca, g, 'show_edges', true, 'show_colorbar', false);

%%
% Create the Fourier transform.

g = grasp_eigendecomposition(g);

%%
% Create the translation operators.

[~, TG] = grasp_translation(g);
T1 = grasp_generalized_translation(g, 1);
GS = g.A;

%% Example graph signal: delta
% Translate a delta signal to study the impulse response.

d10 = grasp_delta(g, 10);
grasp_show_graph(gca, g, 'node_values', abs(d10), 'value_scale', [0 2]);
title('|\delta_{10}|');
snapnow;
grasp_show_graph(gca, g, 'node_values', abs(TG * d10), 'value_scale', [0 1]);
title('|T_g \delta_{10}|');
snapnow;
grasp_show_graph(gca, g, 'node_values', abs(T1 * d10), 'value_scale', [0 max(abs(T1 * d10))], 'highlight_nodes', 1);
title('|T_1 \delta_{10}|');
snapnow;
grasp_show_graph(gca, g, 'node_values', abs(GS * d10), 'value_scale', [0 1]);
title('|A \delta_{10}|');
snapnow;

%% 
% Iterate the graph translation and the graph shift.
% We observe the saturation of the color scale for the graph shift due to
% the non-isometry of this operator.

kmax = 100;
translated_gt(100, kmax + 1) = 0;
translated_gs(100, kmax + 1) = 0;
for k = 0:kmax
    translated_gt(:, k + 1) = TG ^ k * d10;
    translated_gs(:, k + 1) = GS ^ k * d10;
end
mod_options = struct('value_scale', [0 1]);
ang_options = struct('value_scale', [-pi pi], 'color_map', 'hsv');
gs_options = struct('value_scale', [0 1]);
titles_mod = arrayfun(@(k) ['$|T_g^{' int2str(k) '} d_{10}|$'], 0:kmax, 'UniformOutput', false);
titles_ang = arrayfun(@(k) ['$\angle(T_g^{' int2str(k) '} d_{10})$'], 0:kmax, 'UniformOutput', false);
titles_gs = arrayfun(@(k) ['$|A^{' int2str(k) '} d_{10}|$'], 0:kmax, 'UniformOutput', false);
grasp_generate_gif(gcf, 'html/tg_d10_mod_anim.gif', g, abs(translated_gt), titles_mod, 24, mod_options);
grasp_generate_gif(gcf, 'html/tg_d10_ang_anim.gif', g, angle(translated_gt), titles_ang, 24, ang_options);
grasp_generate_gif(gcf, 'html/gs_d10_anim.gif', g, abs(translated_gs), titles_gs, 24, gs_options);
close;

%%
%
% <html>
% <img vspace="5" hspace="5" src="tg_d10_mod_anim.gif" />
% <img vspace="5" hspace="5" src="tg_d10_ang_anim.gif" />
% </html>
%
%%
%
% <html>
% <img vspace="5" hspace="5" src="gs_d10_anim.gif" />
% </html>

%%
% Normalized output for the graph shift

for k = 0:kmax
    translated_gs(:, k + 1) = translated_gs(:, k + 1) / norm(translated_gs(:, k + 1));
end
titles_gs = arrayfun(@(k) ['$|A^{' int2str(k) '} d_{10}| / ||A^{' int2str(k) '} d_{10}||_2$'], 0:kmax, 'UniformOutput', false);
cla(gca);
grasp_generate_gif(gcf, 'html/gs_d10_norm_anim.gif', g, abs(translated_gs), titles_gs, 24, gs_options);
close;

%%
%
% <html>
% <img vspace="5" hspace="5" src="gs_d10_norm_anim.gif" />
% </html>

%% Example graph signal: heat kernel
% As defined in [Shuman, Ricaud, Vandergheynst 2015], $X$ is a heat kernel
% defined as $\widehat{X}(l)=Cexp(-10\lambda_l)$, with $\|X\|_2 = 1$.

cla(gca);
X = grasp_heat_kernel(g, 10);
X = X / norm(X);
grasp_show_graph(gca, g, 'node_values', abs(X));
title('|X|');

%%
% Translate X.

grasp_show_graph(gca, g, 'node_values', abs(TG * X), 'value_scale', [0 max(abs(TG * X))]);
title('|T_g X|');
snapnow;

grasp_show_graph(gca, g, 'node_values', abs(T1 * X), 'value_scale', [0 max(abs(T1 * X))], 'highlight_nodes', 1);
title('|T_1 X|');
snapnow;

grasp_show_graph(gca, g, 'node_values', abs(GS * X), 'value_scale', [0 max(abs(GS * X))]);
title('|A X|');
snapnow;

%% 
% Iterate the graph translation and the graph shift on the signal X.
% The output of the graph shift is normalised.

kmax = 100;
translated_gt(100, kmax + 1) = 0;
translated_gs(100, kmax + 1) = 0;
for k = 0:kmax
    translated_gt(:, k + 1) = TG ^ k * X;
    tmp = GS ^ k * X;
    translated_gs(:, k + 1) = tmp / norm(tmp);
end
mod_options = struct('value_scale', [0 max(abs(translated_gt(:)))]);
ang_options = struct('value_scale', [-pi pi], 'color_map', 'hsv');
gs_options = struct('value_scale', [0 max(abs(translated_gs(:)))]);
titles_mod = arrayfun(@(k) ['$|T_g^{' int2str(k) '} X|$'], 0:kmax, 'UniformOutput', false);
titles_ang = arrayfun(@(k) ['$\angle(T_g^{' int2str(k) '} X)$'], 0:kmax, 'UniformOutput', false);
titles_gs = arrayfun(@(k) ['$|A^{' int2str(k) '} X| / ||A^{' int2str(k) '} X||_2$'], 0:kmax, 'UniformOutput', false);
grasp_generate_gif(gcf, 'html/tg_X_mod_anim.gif', g, abs(translated_gt), titles_mod, 24, mod_options);
grasp_generate_gif(gcf, 'html/tg_X_ang_anim.gif', g, angle(translated_gt), titles_ang, 24, ang_options);
grasp_generate_gif(gcf, 'html/gs_X_norm_anim.gif', g, abs(translated_gs), titles_gs, 24, gs_options);
close;

%%
%
% <html>
% <img vspace="5" hspace="5" src="tg_X_mod_anim.gif" />
% <img vspace="5" hspace="5" src="tg_X_ang_anim.gif" />
% </html>
%
%%
%
% <html>
% <img vspace="5" hspace="5" src="gs_X_norm_anim.gif" />
% </html>
