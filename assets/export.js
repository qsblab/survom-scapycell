window.dash_clientside = Object.assign({}, window.dash_clientside, {
  client: {
    download_cytoscape_png: function(n_clicks) {
      const el = document.getElementById("grn-network-cytoscape");
      if (!el) return "";

      html2canvas(el).then(canvas => {
        const link = document.createElement("a");
        link.download = "grn_network.png";
        link.href = canvas.toDataURL();
        link.click();
      });

      return "";
    }
  }
});
