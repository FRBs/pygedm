import dash
import dash_labs as dl
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

import plotly.express as px
import plotly.graph_objects as pgo

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
import xarray as xr
import h5py

import pygedm

# Load skymap data
dskymap = h5py.File('data/skymap.h5', mode='r')
skymap_dist = dskymap['dist']

skymap_data_ne = xr.DataArray(np.log(dskymap['ne2001/ddm'][:] + 1), dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': dskymap['gl'], 'gb': dskymap['gb']},
                              attrs={'units': 'DM pc/cm3'})
skymap_data_ymw = xr.DataArray(dskymap['ymw16/ddm'], dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': dskymap['gl'], 'gb': dskymap['gb']},
                              attrs={'units': 'DM pc/cm3'})

# APP SETUP
app = dash.Dash(__name__, plugins=[dl.plugins.FlexibleCallbacks()])
app.title = "PyGEDM"
server = app.server 

theme_name = "minty"
css_url = f"https://bootswatch.com/4/{theme_name}/bootstrap.css"

tpl = dl.templates.DbcSidebarTabs(
    app,
    tab_roles=["Plot", "Output", "Skymap"],
    title="PyGEDM: Galactic Electron Density Models",
    theme=css_url, figure_template=True,  sidebar_columns=3,
)

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
        tpl.div_output(role="Output"),
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
        xl = 'Dispersion measure'
        yl = 'Distance'
    else:
        f = pygedm.dist_to_dm
        units = 1.0 * u.kpc
        dmord = float(dmord) * units
        #print(f, dmord)
        yt = 'DM (pc / cm3)'
        xt = 'Distance (kpc)'
        xl = 'Distance'
        yl = 'Dispersion measure'

    if coords == "Galactic (gl, gb)":
        gl = Angle(x0, unit='degree')
        gb = Angle(x1, unit='degree')
        sc = SkyCoord(gl, gb, frame='galactic')
    else:
        ra = Angle(x0, unit='hourangle')
        dec = Angle(x1, unit='degree')
        sc = SkyCoord(ra, dec)
        
    print(sc.galactic.l, sc.galactic.b, dmord, f)
    dout     = f(sc.galactic.l, sc.galactic.b, dmord, method=model)
    dout_ne  = f(sc.galactic.l, sc.galactic.b, dmord, method='ne2001')
    dout_ymw = f(sc.galactic.l, sc.galactic.b, dmord, method='ymw16')

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

    # SKYMAP
    if model == 'NE2001':
        skymap_data = skymap_data_ne
    else:
        skymap_data = skymap_data_ymw
    
    print(skymap_data.shape)
    skymap = px.imshow(animation_frame=0, img=skymap_data,
                        origin='upper', binary_string=True, binary_format='jpg'   )
    skymap["layout"].pop("updatemenus")

    skymap.update_layout(xaxis_title='Galactic longitude (deg)',
                         yaxis_title='Galactic latitude (Deg)',
                         title=f'{model}: All-sky DM map'
                         )
    skymap.add_shape(
        type='rect',
        x0=sc.galactic.l.value-2, x1=sc.galactic.l.value+2, y0=sc.galactic.b.value-2, y1=sc.galactic.b.value+2,
        xref='x', yref='y',
        line_color='cyan'
    )

    ## TEXT OUTPUT
    hdr = html.Div([html.H2(method),
                    html.H4(f"Sky coordinates: {sc}")
          ])

    table_header = [
        html.Thead(html.Tr([html.Th(""), html.Th("YMW16"), html.Th("NE2001")]))
    ]
    row1 = html.Tr([html.Th(f"{xl}"), html.Td(f"{dmord}"), html.Td(f"{dmord}")])
    row2 = html.Tr([html.Th(f"{yl}"), html.Td(f"{dout_ymw[0]:2.4f}"), html.Td(f"{dout_ne[0]:2.4f}")])
    row3 = html.Tr([html.Th("Scattering timescale"), html.Td(f"{dout_ymw[1]:2.4e}"), html.Td(f"{dout_ne[1]:2.4e}")])

    table_body = [html.Tbody([row1, row2, row3])]

    gedm_out = html.Div([hdr, dbc.Table(table_header + table_body, bordered=True)])

    return fig, gedm_out, skymap


app.layout = dbc.Container(fluid=True, children=tpl.children)


if __name__ == "__main__":
    app.run_server(host='0.0.0.0', debug=False)