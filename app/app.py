import dash
import dash_labs as dl
import numpy as np
import dash_core_components as dcc
import plotly.express as px
import plotly.graph_objects as pgo
import pygedm
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
import dash_html_components as html
import xarray as xr

app = dash.Dash(__name__, plugins=[dl.plugins.FlexibleCallbacks()])

theme_name = "minty"
css_url = f"https://bootswatch.com/4/{theme_name}/bootstrap.css"

tpl = dl.templates.DbcSidebarTabs(
    ["Plot", "Output", "Skymap"],
    title="PyGEDM: Galactic Electron Density Models",
    theme=css_url, figure_template=True,  sidebar_columns=3
)

skymap_data_ne = np.load('data/datacube_ddm_ne2001.npy')
skymap_data_ymw = np.load('data/datacube_ddm_ymw16.npy')
skymap_gl = np.load('data/datadim_gl.npy')
skymap_gb = np.load('data/datadim_gb.npy')
skymap_dist = np.round(np.load('data/datadim_dist.npy'), decimals=1)

skymap_data_ne = xr.DataArray(skymap_data_ne, dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': skymap_gl, 'gb': skymap_gb},
                              attrs={'units': 'DM pc/cm3'})
skymap_data_ymw = xr.DataArray(skymap_data_ymw, dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': skymap_gl, 'gb': skymap_gb},
                              attrs={'units': 'DM pc/cm3'})

dist = np.arange(1, 30, 1)
@app.callback(
    args=dict(
        model=tpl.dropdown_input(["NE2001", "YMW16"], label="Model"),
        method=tpl.dropdown_input(["DM (pc/cm3) to Distance", "Distance (kpc) to DM"], label="Method", kind=dl.State),
        dmord=tpl.textbox_input(10, label=None, kind=dl.State),
        coords=tpl.dropdown_input(["Galactic (gl, gb)", "Celestial (RA, DEC)"], label="Coordinates", kind=dl.State),
        x0=tpl.textbox_input("00:00:00.00", label=None, kind=dl.State),
        x1=tpl.textbox_input("00:00:00.00", label=None, kind=dl.State),
        go=tpl.button_input("Calculate", label=None),
        tab=tpl.tab_input(),
    ),
    output=[
        tpl.graph_output(role="Plot"),
        tpl.markdown_output(role="Output"),
        tpl.graph_output(role="Skymap")
    ],
    template=tpl,
)
def callback(model, method, dmord, coords, x0, x1, go, tab):
    print(f"{tab}")
    if method == "DM (pc/cm3) to Distance":
        f = pygedm.dm_to_dist
        units = 1.0 * u.pc / (u.cm**3)
        dmord = float(dmord) * units
        xt = 'DM (pc / cm3)'
        yt = 'Distance (kpc)'
    else:
        f = pygedm.dist_to_dm
        units = 1.0 * u.kpc
        dmord = float(dmord) * units
        print(f, dmord)
        yt = 'DM (pc / cm3)'
        xt = 'Distance (kpc)'

    if coords == "Galactic (gl, gb)":
        gl = Angle(x0, unit='degree')
        gb = Angle(x1, unit='degree')
        sc = SkyCoord(gl, gb, frame='galactic')
        dout = f(gl, gb, dmord, method=model)
    else:
        ra = Angle(x0, unit='hourangle')
        dec = Angle(x1, unit='degree')
        sc = SkyCoord(ra, dec)
        dout = f(sc.galactic.l, sc.galactic.b, dmord, method=model)
        print(sc.galactic.l, sc.galactic.b, dmord, f)

    # Make plots
    D = np.linspace(0.1, dmord.value)
    y_ne21 = np.zeros_like(D)
    y_ymw  = np.zeros_like(D)

    for ii, d in enumerate(D):
        d_ne21 = f(sc.galactic.l, sc.galactic.b, d * units, method='ne2001')
        d_ymw  = f(sc.galactic.l, sc.galactic.b, d * units, method='ymw16')
        if method == "DM (pc/cm3) to Distance":
            y_ne21[ii] = d_ne21[0].to('kpc').value
            y_ymw[ii]  = d_ymw[0].to('kpc').value
        else:
            y_ne21[ii] = d_ne21[0].value
            y_ymw[ii]  = d_ymw[0].value

    #print(d, y)
    fig = pgo.Figure()
    fig.add_trace(pgo.Scatter(x=D, y=y_ne21, mode='lines', name='NE2001'))
    fig.add_trace(pgo.Scatter(x=D, y=y_ymw , mode='lines', name='YMW16'))
    fig.update_layout(xaxis_title=xt, yaxis_title=yt,
                      title=f'ICRS: {sc.icrs.ra:2.2f}  {sc.icrs.dec:2.2f} \t Galactic: {sc.galactic.l:2.2f} {sc.galactic.b:2.2f}')

    if model == 'NE2001':
        skymap_data = skymap_data_ne
    else:
        skymap_data = skymap_data_ymw
    skymap = px.imshow(animation_frame=0, img=skymap_data, x=skymap_gl, y=skymap_gb)
    skymap["layout"].pop("updatemenus")
    skymap.update_layout(xaxis_title='Galactic longitude (deg)',
                         yaxis_title='Galactic latitude (Deg)',
                         title=f'{model}: All-sky DM map')
    skymap.add_shape(
        type='rect',
        x0=gl.value-2, x1=gl.value+2, y0=gb.value-2, y1=gb.value+2,
        xref='x', yref='y',
        line_color='cyan'
    )
    print(skymap)

    gedm_out = f"""# {method} \n
### Input

* **Model**: {model}
* **Sky coordinates:** {sc}
* **DM:** {dmord}

### Output
* **Distance:** {dout[0]:2.4f}
* **Scattering timescale:** {dout[1]:2.4e}
"""

    return fig, gedm_out, skymap


app.layout = tpl.layout(app)


if __name__ == "__main__":
    app.run_server(debug=True)