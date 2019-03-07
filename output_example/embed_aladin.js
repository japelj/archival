function aladin1(nam,ra,dec,fov1,epo){
var aladin = A.aladin('#aladin-lite-div1', {survey: "P/DSS2/red", fov:fov1});
aladin.gotoRaDec(ra, dec);
var marker1=A.marker(ra, dec, {popupTitle: nam});
var markerLayer = A.catalog({sourceSize:18});
aladin.addCatalog(markerLayer);
markerLayer.addSources([marker1]);
aladin.addCatalog(A.catalogFromSimbad({ra: ra, dec: dec}, 0.1, {shape: 'plus', color : '#5d5', onClick: 'showPopup', sourceSize: 18}));
aladin.addCatalog(A.catalogFromVizieR('V/147/sdss12',{ra: ra, dec: dec}, 0.1, {shape: 'square', sourceSize: 8, color: 'red', onClick: 'showPopup'}));
aladin.addCatalog(A.catalogFromVizieR('II/349/ps1',{ra: ra, dec: dec}, 0.1, {shape: 'circle', sourceSize: 8, color: 'blue', onClick: 'showPopup'}));
if (epo > 0.0){
aladin.addCatalog(A.catalogFromSkyBot(ra, dec, 0.1, epo, {shape: 'triangle', sourceSize: 8, color: 'yellow', onClick: 'showPopup'}));
}
}

function aladin2(nam,ra,dec,fov1,epo){
var aladin = A.aladin('#aladin-lite-div2', {survey: "P/PanSTARRS/DR1/color-z-zg-g", fov:fov1});
aladin.gotoRaDec(ra, dec);
var marker1=A.marker(ra, dec, {popupTitle: nam});
var markerLayer = A.catalog({sourceSize:18});
aladin.addCatalog(markerLayer);
markerLayer.addSources([marker1]);
aladin.addCatalog(A.catalogFromSimbad({ra: ra, dec: dec}, 0.1, {shape: 'plus', color : '#5d5', onClick: 'showPopup', sourceSize: 18}));
aladin.addCatalog(A.catalogFromVizieR('V/147/sdss12',{ra: ra, dec: dec}, 0.1, {shape: 'square', sourceSize: 8, color: 'red', onClick: 'showPopup'}));
aladin.addCatalog(A.catalogFromVizieR('II/349/ps1',{ra: ra, dec: dec}, 0.1, {shape: 'circle', sourceSize: 8, color: 'blue', onClick: 'showPopup'}));
if (epo > 0.0){
aladin.addCatalog(A.catalogFromSkyBot(ra, dec, 0.1, epo, {shape: 'triangle', sourceSize: 8, color: 'yellow', onClick: 'showPopup'}));
}
}

function aladin3(nam,ra,dec,fov1,epo){
var aladin = A.aladin('#aladin-lite-div3', {survey: "P/PanSTARRS/DR1/color-z-zg-g", fov:fov1});
aladin.gotoRaDec(ra, dec);
var marker1=A.marker(ra, dec, {popupTitle: nam});
var markerLayer = A.catalog({sourceSize:18});
aladin.addCatalog(markerLayer);
markerLayer.addSources([marker1]);
aladin.addCatalog(A.catalogFromSimbad({ra: ra, dec: dec}, 0.1, {shape: 'plus', color : '#5d5', onClick: 'showPopup', sourceSize: 18}));
aladin.addCatalog(A.catalogFromVizieR('V/147/sdss12',{ra: ra, dec: dec}, 0.1, {shape: 'square', sourceSize: 8, color: 'red', onClick: 'showPopup'}));
aladin.addCatalog(A.catalogFromVizieR('II/349/ps1',{ra: ra, dec: dec}, 0.1, {shape: 'circle', sourceSize: 8, color: 'blue', onClick: 'showPopup'}));
if (epo > 0.0){
aladin.addCatalog(A.catalogFromSkyBot(ra, dec, 0.1, epo, {shape: 'triangle', sourceSize: 8, color: 'yellow', onClick: 'showPopup'}));
}
}

function aladin4(nam,ra,dec,fov1,epo){
var aladin = A.aladin('#aladin-lite-div4', {survey: "P/PanSTARRS/DR1/color-z-zg-g", fov:fov1});
aladin.gotoRaDec(ra, dec);
var marker1=A.marker(ra, dec, {popupTitle: nam});
var markerLayer = A.catalog({sourceSize:18});
aladin.addCatalog(markerLayer);
markerLayer.addSources([marker1]);
aladin.addCatalog(A.catalogFromSimbad({ra: ra, dec: dec}, 0.1, {shape: 'plus', color : '#5d5', onClick: 'showPopup', sourceSize: 18}));
aladin.addCatalog(A.catalogFromVizieR('V/147/sdss12',{ra: ra, dec: dec}, 0.1, {shape: 'square', sourceSize: 8, color: 'red', onClick: 'showPopup'}));
aladin.addCatalog(A.catalogFromVizieR('II/349/ps1',{ra: ra, dec: dec}, 0.1, {shape: 'circle', sourceSize: 8, color: 'blue', onClick: 'showPopup'}));
if (epo > 0.0){
aladin.addCatalog(A.catalogFromSkyBot(ra, dec, 0.1, epo, {shape: 'triangle', sourceSize: 8, color: 'yellow', onClick: 'showPopup'}));
}
}

function aladin5(nam,ra,dec,fov1,epo){
var aladin = A.aladin('#aladin-lite-div5', {survey: "P/PanSTARRS/DR1/color-z-zg-g", fov:fov1});
aladin.gotoRaDec(ra, dec);
var marker1=A.marker(ra, dec, {popupTitle: nam});
var markerLayer = A.catalog({sourceSize:18});
aladin.addCatalog(markerLayer);
markerLayer.addSources([marker1]);
aladin.addCatalog(A.catalogFromSimbad({ra: ra, dec: dec}, 0.1, {shape: 'plus', color : '#5d5', onClick: 'showPopup', sourceSize: 18}));
aladin.addCatalog(A.catalogFromVizieR('V/147/sdss12',{ra: ra, dec: dec}, 0.1, {shape: 'square', sourceSize: 8, color: 'red', onClick: 'showPopup'}));
aladin.addCatalog(A.catalogFromVizieR('II/349/ps1',{ra: ra, dec: dec}, 0.1, {shape: 'circle', sourceSize: 8, color: 'blue', onClick: 'showPopup'}));
if (epo > 0.0){
aladin.addCatalog(A.catalogFromSkyBot(ra, dec, 0.1, epo, {shape: 'triangle', sourceSize: 8, color: 'yellow', onClick: 'showPopup'}));
}
}

