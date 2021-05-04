import dash
import dash_labs as dl
import numpy as np
import dash_core_components as dcc
import plotly.express as px
import plotly.graph_objects as pgo
import pygedm
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u

app = dash.Dash(__name__, plugins=[dl.plugins.FlexibleCallbacks()])

# css_url = "https://bootswatch.com/4/cerulean/bootstrap.css"
#css_url = "https://bootswatch.com/4/cosmo/bootstrap.css"
# css_url = "https://bootswatch.com/4/cyborg/bootstrap.css"
# css_url = "https://bootswatch.com/4/darkly/bootstrap.css"
# css_url = "https://bootswatch.com/4/flatly/bootstrap.css"
# css_url = "https://bootswatch.com/4/journal/bootstrap.css"
# css_url = "https://bootswatch.com/4/litera/bootstrap.css"
# css_url = "https://bootswatch.com/4/lumen/bootstrap.css"
# css_url = "https://bootswatch.com/4/lux/bootstrap.css"
# css_url = "https://bootswatch.com/4/materia/bootstrap.css"
css_url = "https://bootswatch.com/4/minty/bootstrap.css"
# css_url = "https://bootswatch.com/4/pulse/bootstrap.css"
# css_url = "https://bootswatch.com/4/sandstone/bootstrap.css"
# css_url = "https://bootswatch.com/4/simplex/bootstrap.css"
# css_url = "https://bootswatch.com/4/sketchy/bootstrap.css"
# css_url = "https://bootswatch.com/4/slate/bootstrap.css"
# css_url = "https://bootswatch.com/4/solar/bootstrap.css"
# css_url = "https://bootswatch.com/4/spacelab/bootstrap.css"
# css_url = "https://bootswatch.com/4/superhero/bootstrap.css"
# css_url = "https://bootswatch.com/4/united/bootstrap.css"
#css_url = "https://bootswatch.com/4/yeti/bootstrap.css"

tpl = dl.templates.DbcSidebarTabs(
    ["Output", "Galaxy view"],
    title="PyGEDM: Galactic Electron Density Models",
    theme=css_url, figure_template=True,  sidebar_columns=3
)

dist = np.arange(1, 30, 1)
@app.callback(
    args=dict(
        model=tpl.dropdown_input(["NE2001", "YMW16"], label="Model", kind=dl.State),
        method=tpl.dropdown_input(["DM to Distance", "Distance to DM"], label="Method", kind=dl.State),
        dmord=tpl.textbox_input(10, label=None, kind=dl.State),
        coords=tpl.dropdown_input(["Galactic", "Celestial"], label="Coordinates", kind=dl.State),
        x0=tpl.textbox_input("00:00:00.00", label="RA", kind=dl.State),
        x1=tpl.textbox_input("00:00:00.00", label="DEC", kind=dl.State),
        go=tpl.button_input("Calculate", label=None),
        tab=tpl.tab_input(),
        #distance=tpl.slider_input(
        #    dist[0], dist[-1], step=1, value=dist[0], label="Distance"
        #),
        #logs=tpl.checklist_input(
        #    ["DM", "log10(DM)"], value="DM", label="Image Scale"
        #),
    ),
    output=[
        tpl.graph_output(role="Output"),
        tpl.graph_output(role="Galaxy view"),
    ],
    template=tpl,
)
def callback(model, method, dmord, coords, x0, x1, go, tab):
    print(f"{tab}")
    if method == "DM to Distance":
        f = pygedm.dm_to_dist
        dmord = float(dmord) * u.pc / (u.cm**3)
        xt = 'DM (pc / cm3)'
        yt = 'Distance (kpc)'
    else:
        f = pygedm.dist_to_dm
        dmord = float(dmord) * u.kpc
        yt = 'DM (pc / cm3)'
        xt = 'Distance (kpc)'

    if coords == "Galactic":
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
        d_ne21 = f(sc.galactic.l, sc.galactic.b, d, method='ne2001')
        d_ymw  = f(sc.galactic.l, sc.galactic.b, d, method='ymw16')
        print(d_ne21, d_ymw)
        if method == "DM to Distance":
            y_ne21[ii]  = d_ne21[0].to('kpc').value
            y_ymw[ii]   = d_ymw[0].to('kpc').value
        else:
            y_ne21[ii]  = d_ne21[0].value
            y_ymw[ii] = d_ymw[0].value

    #print(d, y)
    fig = pgo.Figure()
    fig.add_trace(pgo.Scatter(x=D, y=y_ne21, mode='lines', name='NE2001'))
    fig.add_trace(pgo.Scatter(x=D, y=y_ymw , mode='lines', name='YMW16'))
    fig.update_layout(xaxis_title=xt, yaxis_title=yt, title=f'{sc}')
    galview = pgo.Figure()

    return (fig, galview)


app.layout = tpl.layout(app)


if __name__ == "__main__":
    app.run_server(debug=True)