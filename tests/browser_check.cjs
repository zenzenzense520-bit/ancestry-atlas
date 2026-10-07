const {chromium}=require('playwright');
const path=require('node:path'),fs=require('node:fs/promises'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..');
(async()=>{
 const expected={};for(const scope of ['global','east_asian'])expected[scope]=JSON.parse(await fs.readFile(path.join(root,`results/${scope}_pca.json`),'utf8'));
 const browser=await chromium.launch({headless:true,executablePath:process.env.ATLAS_CHROMIUM||path.resolve(root,'../browser-local/chromium'),args:['--no-sandbox','--disable-dev-shm-usage','--disable-gpu']});
 const report=[];await fs.mkdir(path.join(root,'qa'),{recursive:true});
 for(const width of [360,390,768,1440]){
  const context=await browser.newContext({viewport:{width,height:900},acceptDownloads:true});const page=await context.newPage(),errors=[],requests=[];
  page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url());});
  await page.goto('file://'+path.join(root,'report/Sample01_reference_report_v05.html'));await page.evaluate(()=>document.fonts.ready);
  const dimensions=await page.evaluate(()=>({root:document.documentElement.scrollWidth,viewport:innerWidth}));assert.ok(dimensions.root<=width+1,JSON.stringify(dimensions));
  const ids=await page.locator('[id]').evaluateAll(e=>e.map(x=>x.id));assert.equal(new Set(ids).size,ids.length);
  for(const scope of ['global','east_asian']){
   const module=page.locator(`[data-scope="${scope}"]`),p=expected[scope];
   const counts=await module.locator('svg [data-group]').evaluateAll(groups=>groups.map(g=>g.querySelectorAll('use').length));
   assert.equal(counts.reduce((a,b)=>a+b,0),p.reference.reference_n*4);
   await module.locator('[data-view="3"]').click();assert.equal(await module.locator('[data-chart="3"]').isVisible(),true);
   assert.match(await module.locator('.view-caption').innerText(),new RegExp((p.variance_explained[2]*100).toFixed(2).replace('.','\\.')));
   await module.locator('[data-view="3"]').focus();await page.keyboard.press('ArrowLeft');assert.equal(await module.locator('[data-chart="2"]').isVisible(),true);
   const group=Object.keys(p.reference.group_counts)[0];const toggle=module.locator(`[data-group-toggle="${group}"]`);
   await toggle.click();assert.equal(await module.locator(`svg [data-group="${group}"].hidden-group`).count(),4);
   await toggle.focus();await page.keyboard.press('Space');assert.equal(await module.locator('svg .hidden-group').count(),0);
  }
  await page.locator('[data-scope="global"] [data-view="3"]').click();assert.equal(await page.locator('[data-scope="east_asian"] [data-chart="2"]').isVisible(),true);
  await page.locator('[data-scope="global"] [data-view="2"]').click();
  const promise=page.waitForEvent('download');await page.locator('#export-summary').click();const download=await promise;
  const payload=JSON.parse(await fs.readFile(await download.path(),'utf8'));
  assert.equal(payload.global.used_after_ld,expected.global.used_after_ld);assert.equal(payload.east_asian.nearest_reference_center,expected.east_asian.nearest_reference_center);
  assert.ok(!('genotypes' in payload));
  await page.locator('.chromosome-detail summary').click();assert.equal(await page.locator('.chromosome-table').isVisible(),true);await page.locator('.chromosome-detail summary').click();
  await page.evaluate(()=>scrollTo({top:0,behavior:'instant'}));await page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))));
  assert.equal(await page.locator('nav a[aria-current]').getAttribute('href'),'#finding');
  await page.screenshot({path:path.join(root,`qa/report-${width}.png`),fullPage:true,animations:'disabled'});
  if(width===1440||width===390){
   await page.screenshot({path:path.join(root,`qa/hero-${width}.png`),animations:'disabled'});
   await page.evaluate(()=>scrollTo({top:document.getElementById('east').offsetTop-110,behavior:'instant'}));
   await page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))));
   await page.screenshot({path:path.join(root,`qa/east-${width}.png`),animations:'disabled'});
   assert.equal(await page.locator('nav a[aria-current]').getAttribute('href'),'#east');
  }
  assert.deepEqual(errors,[]);assert.deepEqual(requests,[]);
  report.push({width,noOverflow:true,uniqueIds:ids.length,globalPoints:expected.global.reference.reference_n,eastAsianPoints:expected.east_asian.reference.reference_n,independentViews:true,keyboard:true,jsonExport:true,noExternalRequests:true,noErrors:true});await context.close();
 }
 const context=await browser.newContext({viewport:{width:390,height:844},reducedMotion:'reduce'});const page=await context.newPage();await page.goto('file://'+path.join(root,'report/Sample01_reference_report_v05.html'));
 assert.equal(await page.locator('h1').evaluate(e=>getComputedStyle(e).animationName),'none');await page.locator('[data-scope="global"] [data-group-toggle="EAS"]').click();
 await page.emulateMedia({media:'print'});assert.equal(await page.locator('.topbar').isVisible(),false);
 assert.equal(await page.locator('[data-scope="global"] [data-chart="2"] .desktop-chart').isVisible(),true);
 assert.equal(await page.locator('[data-scope="global"] [data-chart="2"] .mobile-chart').isVisible(),false);
 assert.equal(await page.locator('svg .hidden-group').first().evaluate(e=>getComputedStyle(e).opacity),'1');
 report.push({reducedMotion:true,printMode:true,printRestoresPoints:true});await context.close();await browser.close();
 await fs.writeFile(path.join(root,'qa/browser_checks.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exit(1);});
