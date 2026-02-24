import dash
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, ctx, no_update
import h5py
import numpy as np
import plotly.express as px
import plotly.graph_objects as pgo
import xarray as xr
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

import pygedm

# Load skymap data
dskymap = h5py.File("data/skymap.h5", mode="r")
skymap_dist = dskymap["dist"]

skymap_data_ne = xr.DataArray(
    dskymap["ne2001/ddm"][:, ::2, ::2],
    dims=("distance_kpc", "gb", "gl"),
    coords={
        "distance_kpc": skymap_dist,
        "gl": dskymap["gl"][::2],
        "gb": dskymap["gb"][::2],
    },
    attrs={"units": "DM pc/cm3"},
)
skymap_data_ymw = xr.DataArray(
    dskymap["ymw16/ddm"][:, ::2, ::2],
    dims=("distance_kpc", "gb", "gl"),
    coords={
        "distance_kpc": skymap_dist,
        "gl": dskymap["gl"][::2],
        "gb": dskymap["gb"][::2],
    },
    attrs={"units": "DM pc/cm3"},
)

# APP SETUP
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SPACELAB])
app.title = "PyGEDM"
server = app.server


@app.callback(
    Output("plot-output", "figure"),
    Output("table-output", "children"),
    Output("skymap-output", "figure"),
    Output("notes-output", "children"),
    Output("dm-min-input", "value"),
    Output("dm-max-input", "value"),
    Output("skymap-slider-store", "data"),
    Input("calculate-button", "n_clicks"),
    Input("skymap-apply-button", "n_clicks"),
    Input("skymap-output", "relayoutData"),
    State("model-dropdown", "value"),
    State("colorscale-dropdown", "value"),
    State("dm-min-input", "value"),
    State("dm-max-input", "value"),
    State("skymap-output", "figure"),
    State("skymap-slider-store", "data"),
    State("method-dropdown", "value"),
    State("dmord-input", "value"),
    State("nu-input", "value"),
    State("coords-dropdown", "value"),
    State("x0-input", "value"),
    State("x1-input", "value"),
    prevent_initial_call=False,
)
def callback(n_clicks, skymap_apply_clicks, relayout_data, model, colorscale, dm_min, dm_max, skymap_fig, slider_store, method, dmord, nu, coords, x0, x1):
    """Main callback for all calculations and plots"""

    # Setup some error handling
    coord_error = False
    freq_error = False
    dmord_error = False

    # Check what triggered the callback
    triggered_id = ctx.triggered_id if ctx.triggered else None
    initial_load = (
        (not ctx.triggered)
        or (ctx.triggered[0].get("prop_id") in (".", ""))
        or (triggered_id in (None, "", "."))
        or (
            triggered_id == "skymap-output"
            and relayout_data == {"autosize": True}
            and (n_clicks in (None, 0))
            and (skymap_apply_clicks in (None, 0))
        )
    )

    # Check if this is just a zoom/pan event - if so, don't regenerate the skymap
    is_zoom_or_pan = False
    if triggered_id == 'skymap-output' and relayout_data:
        # Check if it's a zoom/pan (contains axis ranges)
        if any(k in relayout_data for k in ['xaxis.range', 'yaxis.range', 'xaxis.range[0]', 'yaxis.range[0]',
                                              'xaxis.autorange', 'yaxis.autorange']):
            is_zoom_or_pan = True

    # Check if distance slider changed (relayoutData contains frame updates)
    reset_dm_range = False
    current_frame_idx = 0
    slider_store_out = no_update
    default_slider_idx = None
    slider_value_candidate = None
    if relayout_data and not is_zoom_or_pan:
        # Check if slider.value key exists (frame changed)
        if 'slider.value' in relayout_data:
            reset_dm_range = True
            current_frame_idx = relayout_data['slider.value']
            slider_value_candidate = current_frame_idx
        # Also check for the frame argument pattern
        elif any('frame' in str(k).lower() for k in relayout_data.keys()):
            reset_dm_range = True
            # Try to capture frame value if present
            if 'frame' in relayout_data:
                slider_value_candidate = relayout_data['frame']

        # Try other slider keys that plotly may emit
        if slider_value_candidate is None:
            for key in (
                'sliders[0].active',
                'sliders[0].value',
                'slider.active',
                'slider.value',
            ):
                if key in relayout_data:
                    slider_value_candidate = relayout_data[key]
                    break

    # Decide which outputs to update
    update_plot = initial_load or triggered_id == "calculate-button"
    update_skymap = initial_load or triggered_id == "skymap-apply-button" or (
        triggered_id == "skymap-output" and reset_dm_range
    )

    # Capture current slider active index from existing figure (e.g., when Apply is clicked)
    if slider_value_candidate is None and update_skymap and isinstance(skymap_fig, dict):
        try:
            sliders = skymap_fig.get("layout", {}).get("sliders", [])
            if sliders:
                slider_value_candidate = sliders[0].get("active")
        except AttributeError:
            slider_value_candidate = None

    if slider_value_candidate is not None:
        slider_store_out = slider_value_candidate

    try:
        nu = float(nu)
    except ValueError:
        nu = 1.0
        freq_error = True

    if method == "DM (pc/cm3) to Distance":
        f = pygedm.dm_to_dist
        units = 1.0 * u.pc / (u.cm**3)
        try:
            dmord = float(dmord) * units
        except ValueError:
            dmord = 10 * units
            dmord_error = True
        xt = "DM (pc / cm3)"
        yt = "Distance (kpc)"
        xl = "Dispersion measure"
        yl = "Distance"
    else:
        f = pygedm.dist_to_dm
        units = 1.0 * u.kpc
        try:
            dmord = float(dmord) * units
        except ValueError:
            dmord = 10 * units
            dmord_error = True
        # print(f, dmord)
        yt = "DM (pc / cm3)"
        xt = "Distance (kpc)"
        xl = "Distance"
        yl = "Dispersion measure"

    if coords == "Galactic (gl, gb)":
        try:
            gl = Angle(x0, unit="degree")
            gb = Angle(x1, unit="degree")
            sc = SkyCoord(gl, gb, frame="galactic")
        except ValueError:
            sc = SkyCoord(0 * u.deg, 0 * u.deg, frame="galactic")
            coord_error = True
    else:
        try:
            ra = Angle(x0, unit="hourangle")
            dec = Angle(x1, unit="degree")
            sc = SkyCoord(ra, dec)
        except:
            sc = SkyCoord(0 * u.deg, 0 * u.deg, frame="galactic")
            coord_error = True

    print(sc.galactic.l, sc.galactic.b, dmord, f)
    dout = f(sc.galactic.l, sc.galactic.b, dmord, method=model, nu=nu)
    dout_ne = f(sc.galactic.l, sc.galactic.b, dmord, method="ne2001", nu=nu)
    dout_ne25 = f(sc.galactic.l, sc.galactic.b, dmord, method="ne2025", nu=nu)
    dout_ymw = f(sc.galactic.l, sc.galactic.b, dmord, method="ymw16", nu=nu)

    # Make plots
    D = np.linspace(0.1, dmord.value)
    y_ne21 = np.zeros_like(D)
    y_ne25 = np.zeros_like(D)
    y_ymw = np.zeros_like(D)

    for ii, d in enumerate(D):
        d_ne21 = f(sc.galactic.l, sc.galactic.b, d * units, method="ne2001", nu=nu)
        d_ne25 = f(sc.galactic.l, sc.galactic.b, d * units, method="ne2025", nu=nu)
        d_ymw = f(sc.galactic.l, sc.galactic.b, d * units, method="ymw16", nu=nu)
        if method == "DM (pc/cm3) to Distance":
            y_ne21[ii] = d_ne21[0].to("kpc").value
            y_ne25[ii] = d_ne25[0].to("kpc").value
            y_ymw[ii] = d_ymw[0].to("kpc").value
        else:
            y_ne21[ii] = d_ne21[0].value
            y_ne25[ii] = d_ne25[0].value
            y_ymw[ii] = d_ymw[0].value

    # print(d, y)
    fig = pgo.Figure()
    fig.add_trace(pgo.Scatter(x=D, y=y_ne21, mode="lines", name="NE2001"))
    fig.add_trace(pgo.Scatter(x=D, y=y_ne25, mode="lines", name="NE2025"))
    fig.add_trace(pgo.Scatter(x=D, y=y_ymw, mode="lines", name="YMW16"))
    fig.update_layout(
        xaxis_title=xt,
        yaxis_title=yt,
        title=f"ICRS: {sc.icrs.ra:2.2f}  {sc.icrs.dec:2.2f} \t Galactic: {sc.galactic.l:2.2f} {sc.galactic.b:2.2f}",
    )

    # SKYMAP
    # NE2025 uses NE2001 data for now (NE2025 skymap not yet available)
    if model == "NE2001" or model == "NE2025":
        skymap_data = skymap_data_ne
        skymap_model_label = "NE2001" if model == "NE2001" else "NE2001 (NE2025 skymap not yet available)"
    else:
        skymap_data = skymap_data_ymw
        skymap_model_label = model

    print(skymap_data.shape)

    # Determine DM min/max for colorscale
    # Only reset if distance slider changed
    if reset_dm_range:
        # Get the data range from the current frame only
        current_frame_data = skymap_data.values[current_frame_idx, :, :]
        dm_min_val = float(np.nanmin(current_frame_data))
        dm_max_val = float(np.nanmax(current_frame_data))
    else:
        # Keep user's values and don't update the input fields
        def _to_float(val):
            try:
                return float(val)
            except (TypeError, ValueError):
                return None

        dm_min_val = _to_float(dm_min)
        dm_max_val = _to_float(dm_max)

    # Determine effective range for color scale (use active frame if one bound is missing)
    active_idx_for_range = current_frame_idx
    if slider_store is not None:
        try:
            if isinstance(slider_store, (int, np.integer)):
                active_idx_for_range = int(slider_store)
            else:
                dist_vals = np.asarray(skymap_data["distance_kpc"].values, dtype=float)
                active_val = float(slider_store)
                active_idx_for_range = int(np.argmin(np.abs(dist_vals - active_val)))
        except (TypeError, ValueError):
            pass

    if dm_min_val is None or dm_max_val is None:
        frame_data = skymap_data.values[active_idx_for_range, :, :]
        if dm_min_val is None:
            dm_min_val = float(np.nanmin(frame_data))
        if dm_max_val is None:
            dm_max_val = float(np.nanmax(frame_data))

    range_color = [dm_min_val, dm_max_val] if (dm_min_val is not None or dm_max_val is not None) else None

    skymap = px.imshow(
        skymap_data,
        origin="lower",
        animation_frame="distance_kpc",
        color_continuous_scale=colorscale,
        range_color=range_color,
    )
    # Remove animation menus if present
    if "updatemenus" in skymap["layout"]:
        skymap["layout"].pop("updatemenus")

    # Default/preserve slider value
    if "sliders" in skymap["layout"] and skymap["layout"]["sliders"]:
        dist_vals = np.asarray(skymap_data["distance_kpc"].values, dtype=float)

        if slider_store is not None and update_skymap:
            try:
                # If stored value is an index
                if isinstance(slider_store, (int, np.integer)):
                    active_idx = int(slider_store)
                else:
                    # Try to parse stored distance value
                    active_val = float(slider_store)
                    active_idx = int(np.argmin(np.abs(dist_vals - active_val)))
                skymap["layout"]["sliders"][0]["active"] = active_idx
            except (ValueError, TypeError):
                pass
        elif update_skymap:
            default_slider_idx = int(np.argmin(np.abs(dist_vals - 8.5)))
            skymap["layout"]["sliders"][0]["active"] = default_slider_idx
            if slider_store is None and slider_store_out is no_update:
                slider_store_out = default_slider_idx

    if initial_load and slider_store is None and slider_store_out is no_update and default_slider_idx is not None:
        slider_store_out = default_slider_idx

    # Preserve zoom/pan state if it exists in relayout_data
    xaxis_range = None
    yaxis_range = None
    if relayout_data and not reset_dm_range:
        # Extract current axis ranges if they exist
        if 'xaxis.range[0]' in relayout_data and 'xaxis.range[1]' in relayout_data:
            xaxis_range = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
        elif 'xaxis.range' in relayout_data:
            xaxis_range = relayout_data['xaxis.range']

        if 'yaxis.range[0]' in relayout_data and 'yaxis.range[1]' in relayout_data:
            yaxis_range = [relayout_data['yaxis.range[0]'], relayout_data['yaxis.range[1]']]
        elif 'yaxis.range' in relayout_data:
            yaxis_range = relayout_data['yaxis.range']

    skymap.update_layout(
        xaxis_title="Galactic longitude (deg)",
        yaxis_title="Galactic latitude (Deg)",
        title=f"{skymap_model_label}: All-sky DM map",
        height=600,
        coloraxis_colorbar=dict(title="DM (pc/cm³)"),
    )

    # Apply preserved zoom if available
    if xaxis_range is not None:
        skymap.update_xaxes(range=xaxis_range)
    if yaxis_range is not None:
        skymap.update_yaxes(range=yaxis_range)

    l_wrapped = sc.galactic.l.wrap_at(Angle(180, unit="deg")).value

    skymap.add_shape(
        type="rect",
        x0=l_wrapped - 2,
        x1=l_wrapped + 2,
        y0=sc.galactic.b.value - 2,
        y1=sc.galactic.b.value + 2,
        xref="x",
        yref="y",
        line_color="cyan",
    )

    ## TEXT OUTPUT
    hdr = [html.H2(method)]
    if coord_error:
        hdr.append(dbc.Alert("Input coordinates invalid, please check", color="danger"))
    if freq_error:
        hdr.append(
            dbc.Alert("Could not parse frequency input, please check.", color="danger")
        )
    if dmord_error:
        hdr.append(
            dbc.Alert(
                "Could not parse DM/distance input, please check.", color="danger"
            )
        )
    hdr = html.Div(hdr)

    table_header = [
        html.Thead(html.Tr([html.Th(""), html.Th("YMW16"), html.Th("NE2025"), html.Th("NE2001")]))
    ]
    row1 = html.Tr([html.Th(f"{xl}"), html.Td(f"{dmord}"), html.Td(f"{dmord}"), html.Td(f"{dmord}")])
    row2 = html.Tr([
        html.Th(f"{yl}"),
        html.Td(f"{dout_ymw[0]:2.4f}"),
        html.Td(f"{dout_ne25[0]:2.4f}"),
        html.Td(f"{dout_ne[0]:2.4f}"),
    ])
    row3 = html.Tr([
        html.Th(f"Scattering timescale @ {nu} GHz"),
        html.Td(f"{dout_ymw[1]:2.4e}"),
        html.Td(f"{dout_ne25[1]:2.4e}"),
        html.Td(f"{dout_ne[1]:2.4e}"),
    ])

    table_body = [html.Tbody([row1, row2, row3])]

    gedm_out = html.Div([hdr, dbc.Table(table_header + table_body, bordered=True)])

    notes = html.Div([
        html.H2("About"),
        html.P([html.B("Version: "), f"{pygedm.__version__}"]),
        html.P([html.B("Authors: "), "Danny C. Price, Chris Flynn,  Adam Deller"]),
        html.P([
            html.B("PyGEDM Documentation: "),
            html.A(
                "https://pygedm.readthedocs.io", href="https://pygedm.readthedocs.io"
            ),
        ]),
        html.P([
            html.B("Github: "),
            html.A(
                "https://github.com/FRBs/pygedm", href="https://github.com/FRBs/pygedm"
            ),
        ]),
        html.H2("References"),
        html.H4("PyGEDM"),
        html.P("Price, D. C., Flynn, C., and Deller, A."),
        html.P([
            html.A(
                "A comparison of Galactic electron density models using PyGEDM",
                href="https://scixplorer.org/abs/2021PASA...38...38P/abstract",
            )
        ]),
        html.H4("NE2001"),
        html.P("Cordes, J. M., & Lazio, T. J. W. (2002),"),
        html.P([
            html.A(
                "NE2001.I. A New Model for the Galactic Distribution of Free Electrons and its Fluctuations, arXiv e-prints, astro-ph/0207156.",
                href="https://ui.adsabs.harvard.edu/abs/2002astro.ph..7156C/abstract"
            )
        ]),
        html.H4("NE2025"),
        html.P("Ocker, S.K. and Cordes, J.M. (2026),"),
        html.P([
            html.A(
                "NE2025: An Updated Electron Density Model for the Galactic Interstellar Medium, arXiv e-prints, astro-ph/2602.11838.",
                href="https://ui.adsabs.harvard.edu/abs/2026arXiv260211838O/abstract"
            )
        ]),
        html.H4("YMW16"),
        html.P("Yao, J. M., Manchester, R. N., & Wang, N. (2017),"),
        html.P([
            html.A("A New Electron-density Model for Estimation of Pulsar and FRB Distances, ApJ, Volume 888, Issue 2, id.105, Col 835, id.29",
                   href="https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract"
            )
        ])
    ])

    # Handle different update scenarios
    if is_zoom_or_pan:
        # Don't regenerate anything on zoom/pan, just let the figure handle it
        return no_update, no_update, no_update, no_update, no_update, no_update, slider_store_out

    plot_out = fig if update_plot else no_update
    table_out = gedm_out if update_plot else no_update
    notes_out = notes if update_plot else no_update

    skymap_out = skymap if update_skymap else no_update
    dm_min_out = dm_min_val if (update_skymap and reset_dm_range) else no_update
    dm_max_out = dm_max_val if (update_skymap and reset_dm_range) else no_update

    return plot_out, table_out, skymap_out, notes_out, dm_min_out, dm_max_out, slider_store_out


# APP LAYOUT
app.layout = dbc.Container(
    fluid=True,
    children=[
        dcc.Store(id="skymap-slider-store"),
        html.H1("PyGEDM: Galactic Electron Density Models", className="my-4"),
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Settings", className="card-title"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Method"),
                                dcc.Dropdown(
                                    id="method-dropdown",
                                    options=[
                                        {"label": "DM (pc/cm3) to Distance", "value": "DM (pc/cm3) to Distance"},
                                        {"label": "Distance (kpc) to DM", "value": "Distance (kpc) to DM"},
                                    ],
                                    value="DM (pc/cm3) to Distance",
                                    searchable=False,
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("DM / Distance Value"),
                                dbc.Input(
                                    id="dmord-input",
                                    type="number",
                                    value=10,
                                    placeholder="Enter value",
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Frequency (GHz)"),
                                dbc.Input(
                                    id="nu-input",
                                    type="number",
                                    value=1.0,
                                    placeholder="1.0",
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Coordinates"),
                                dcc.Dropdown(
                                    id="coords-dropdown",
                                    options=[
                                        {"label": "Galactic (gl, gb)", "value": "Galactic (gl, gb)"},
                                        {"label": "Celestial (RA, DEC)", "value": "Celestial (RA, DEC)"},
                                    ],
                                    value="Galactic (gl, gb)",
                                    searchable=False,
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("X (gl or RA)"),
                                dbc.Input(
                                    id="x0-input",
                                    type="text",
                                    value="00:00:00.00",
                                    placeholder="00:00:00.00",
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Y (gb or DEC)"),
                                dbc.Input(
                                    id="x1-input",
                                    type="text",
                                    value="00:00:00.00",
                                    placeholder="00:00:00.00",
                                ),
                            ], width=12, className="mt-2"),
                        ]),
                        dbc.Row([
                            dbc.Col([
                                dbc.Button(
                                    "Calculate",
                                    id="calculate-button",
                                    color="primary",
                                    className="w-100 mt-3",
                                    n_clicks=0,
                                ),
                            ], width=12),
                        ]),
                    ])
                ])
            ], width=3),
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Output", children=[
                        dcc.Loading(
                            id="loading",
                            type="default",
                            children=[
                                dcc.Graph(id="plot-output"),
                                html.Div(id="table-output"),
                            ],
                        ),
                    ]),
                    dbc.Tab(label="Skymap", children=[
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Model"),
                                dcc.Dropdown(
                                    id="model-dropdown",
                                    options=[
                                        {"label": "NE2001", "value": "NE2001"},
                                        {"label": "NE2025", "value": "NE2025"},
                                        {"label": "YMW16", "value": "YMW16"},
                                    ],
                                    value="YMW16",
                                    searchable=False,
                                    style={"maxWidth": "200px"},
                                ),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("Color Scale"),
                                dcc.Dropdown(
                                    id="colorscale-dropdown",
                                    options=[
                                        {"label": "Viridis", "value": "Viridis"},
                                        {"label": "Plasma", "value": "Plasma"},
                                        {"label": "Inferno", "value": "Inferno"},
                                        {"label": "Magma", "value": "Magma"},
                                        {"label": "Cividis", "value": "Cividis"},
                                        {"label": "Hot", "value": "Hot"},
                                        {"label": "Blues", "value": "Blues"},
                                        {"label": "Greys", "value": "Greys"},
                                    ],
                                    value="Viridis",
                                    searchable=False,
                                    style={"maxWidth": "200px"},
                                ),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("DM Min"),
                                dbc.Input(
                                    id="dm-min-input",
                                    type="number",
                                    placeholder="Auto",
                                    style={"maxWidth": "150px"},
                                ),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("DM Max"),
                                dbc.Input(
                                    id="dm-max-input",
                                    type="number",
                                    placeholder="Auto",
                                    style={"maxWidth": "150px"},
                                ),
                            ], width=2),
                            dbc.Col([
                                dbc.Label(" "),  # Empty label for alignment
                                dbc.Button(
                                    "Apply",
                                    id="skymap-apply-button",
                                    color="primary",
                                    className="w-100",
                                    n_clicks=0,
                                ),
                            ], width=1),
                        ], style={"marginBottom": "1rem"}, align="center"),
                        dcc.Loading(
                            id="loading-skymap",
                            type="default",
                            children=[
                                dcc.Graph(id="skymap-output"),
                            ],
                        ),
                    ]),
                    dbc.Tab(label="About", children=[
                        html.Div(style={"marginTop": "2rem"}, children=html.Div(id="notes-output")),
                    ]),
                ]),
            ], width=9),
        ]),
        html.Hr(style={"marginTop": "3rem"}),
        html.Footer(
            html.P(
                ["Web app hosted by ", html.A("Data Central", href="https://datacentral.org.au/")],
                style={"textAlign": "center", "color": "#666", "fontSize": "0.9rem", "marginBottom": "1rem"}
            )
        ),
    ],
)


if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=False)
