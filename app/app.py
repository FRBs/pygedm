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

skymap_data_ne = xr.DataArray(dskymap['ne2001/ddm'][:, ::2, ::2], dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': dskymap['gl'][::2], 'gb': dskymap['gb'][::2]},
                              attrs={'units': 'DM pc/cm3'})
skymap_data_ymw = xr.DataArray(dskymap['ymw16/ddm'][:, ::2, ::2], dims=('distance_kpc', 'gb', 'gl'),
                              coords={'distance_kpc': skymap_dist, 'gl': dskymap['gl'][::2], 'gb': dskymap['gb'][::2]},
                              attrs={'units': 'DM pc/cm3'})

# APP SETUP
app = dash.Dash(__name__, plugins=[dl.plugins.FlexibleCallbacks()])
app.title = "PyGEDM"
server = app.server 

theme_name = "minty"
css_url = f"https://bootswatch.com/4/{theme_name}/bootstrap.css"

tpl = dl.templates.DbcSidebarTabs(
    app,
    tab_locations=["Plot", "Output", "Skymap"],
    title="PyGEDM: Galactic Electron Density Models",
    theme=css_url, figure_template=True,  sidebar_columns=3,
)

@app.callback(
    args=dict(
        model=tpl.new_dropdown(["NE2001", "YMW16"], label="Model"),
        method=tpl.new_dropdown(["DM (pc/cm3) to Distance", "Distance (kpc) to DM"], label="Method", kind=dl.State),
        dmord=tpl.new_textbox(10, label=None, kind=dl.State),
        nu=tpl.new_textbox(1.0, label='Frequency (GHz)', kind=dl.State),
        coords=tpl.new_dropdown(["Galactic (gl, gb)", "Celestial (RA, DEC)"], label="Coordinates", kind=dl.State),
        x0=tpl.new_textbox("00:00:00.00", label=None, kind=dl.State),
        x1=tpl.new_textbox("00:00:00.00", label=None, kind=dl.State),
        go=tpl.new_button("Calculate", label=None),
        tab=tpl.tab_input(),
    ),
    output=[
        tpl.new_graph(location="Plot"),
        tpl.new_div(location="Output"),
        tpl.new_graph(location="Skymap")
    ],
    template=tpl,
)
def callback(model, method, dmord, nu, coords, x0, x1, go, tab):
    print(f"{tab}")
    
    # Setup some error handling
    coord_error = False
    freq_error  = False
    dmord_error = False
    
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
        xt = 'DM (pc / cm3)'
        yt = 'Distance (kpc)'
        xl = 'Dispersion measure'
        yl = 'Distance'
    else:
        f = pygedm.dist_to_dm
        units = 1.0 * u.kpc
        try:
            dmord = float(dmord) * units
        except ValueError:
            dmord = 10 * units
            dmord_error = True
        #print(f, dmord)
        yt = 'DM (pc / cm3)'
        xt = 'Distance (kpc)'
        xl = 'Distance'
        yl = 'Dispersion measure'

    if coords == "Galactic (gl, gb)":
        try:
            gl = Angle(x0, unit='degree')
            gb = Angle(x1, unit='degree')
            sc = SkyCoord(gl, gb, frame='galactic')
        except ValueError:
            sc = SkyCoord(0*u.deg, 0*u.deg, frame='galactic')
            coord_error = True
    else:
        try:
            ra = Angle(x0, unit='hourangle')
            dec = Angle(x1, unit='degree')
            sc = SkyCoord(ra, dec)
        except:
            sc = SkyCoord(0*u.deg, 0*u.deg, frame='galactic')
            coord_error = True
            
    print(sc.galactic.l, sc.galactic.b, dmord, f)
    dout     = f(sc.galactic.l, sc.galactic.b, dmord, method=model, nu=nu)
    dout_ne  = f(sc.galactic.l, sc.galactic.b, dmord, method='ne2001', nu=nu)
    dout_ymw = f(sc.galactic.l, sc.galactic.b, dmord, method='ymw16', nu=nu)

    # Make plots
    D = np.linspace(0.1, dmord.value)
    y_ne21 = np.zeros_like(D)
    y_ymw  = np.zeros_like(D)

    for ii, d in enumerate(D):
        d_ne21 = f(sc.galactic.l, sc.galactic.b, d * units, method='ne2001', nu=nu)
        d_ymw  = f(sc.galactic.l, sc.galactic.b, d * units, method='ymw16', nu=nu)
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
    skymap = px.imshow(animation_frame=0, img=skymap_data, origin='upper' )
    skymap["layout"].pop("updatemenus")

    skymap.update_layout(xaxis_title='Galactic longitude (deg)',
                         yaxis_title='Galactic latitude (Deg)',
                         title=f'{model}: All-sky DM map'
                         )
                         
    l_wrapped = sc.galactic.l.wrap_at(Angle(180, unit='deg')).value
    
    skymap.add_shape(
        type='rect',
        x0=l_wrapped-2, x1=l_wrapped+2, y0=sc.galactic.b.value-2, y1=sc.galactic.b.value+2,
        xref='x', yref='y',
        line_color='cyan'
    )

    ## TEXT OUTPUT
    hdr = [html.H2(method), html.H4(f"Sky coordinates: {sc}")]
    if coord_error:
        hdr.append(dbc.Alert('Input coordinates invalid, please check', color='danger'))
    if freq_error:
        hdr.append(dbc.Alert('Could not parse frequency input, please check.', color='danger'))
    if dmord_error:
        hdr.append(dbc.Alert('Could not parse DM/distance input, please check.', color='danger'))        
    hdr = html.Div(hdr)

    table_header = [
        html.Thead(html.Tr([html.Th(""), html.Th("YMW16"), html.Th("NE2001")]))
    ]
    row1 = html.Tr([html.Th(f"{xl}"), html.Td(f"{dmord}"), html.Td(f"{dmord}")])
    row2 = html.Tr([html.Th(f"{yl}"), html.Td(f"{dout_ymw[0]:2.4f}"), html.Td(f"{dout_ne[0]:2.4f}")])
    row3 = html.Tr([html.Th(f"Scattering timescale @ {nu} GHz"), html.Td(f"{dout_ymw[1]:2.4e}"), html.Td(f"{dout_ne[1]:2.4e}")])

    table_body = [html.Tbody([row1, row2, row3])]

    gedm_out = html.Div([hdr, dbc.Table(table_header + table_body, bordered=True)])

    return fig, gedm_out, skymap


app.layout = dbc.Container(fluid=True, children=tpl.children)


if __name__ == "__main__":
    app.run_server(host='0.0.0.0', debug=False)