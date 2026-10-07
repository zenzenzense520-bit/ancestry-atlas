'use strict';
(() => {
  const result = JSON.parse(document.getElementById('sample-data').textContent);
  const status = document.getElementById('io-status');
  const reduce = matchMedia('(prefers-reduced-motion: reduce)');
  document.querySelectorAll('[data-scope]').forEach(module => {
  module.querySelectorAll('[data-view]').forEach(button => button.addEventListener('click', () => {
    const view = button.dataset.view;
    module.querySelectorAll('[data-view]').forEach(b => b.setAttribute('aria-pressed', String(b === button)));
    module.querySelectorAll('[data-chart]').forEach(pane => pane.hidden = pane.dataset.chart !== view);
    const variance = result[module.dataset.scope].variance_explained;
    module.querySelector('.view-caption').textContent = `PC1 ${(variance[0]*100).toFixed(2)}% · PC${view} ${(variance[Number(view)-1]*100).toFixed(2)}% 解释方差`;
  }));
  module.querySelectorAll('[data-view]').forEach(button => button.addEventListener('keydown',event=>{
    if(!['ArrowLeft','ArrowRight'].includes(event.key))return;
    event.preventDefault();
    const buttons=[...module.querySelectorAll('[data-view]')];
    const next=buttons[(buttons.indexOf(button)+1)%buttons.length];next.focus();next.click();
  }));
  module.querySelectorAll('[data-group-toggle]').forEach(button => button.addEventListener('click', () => {
    const hide = button.getAttribute('aria-pressed') === 'true';
    button.setAttribute('aria-pressed', String(!hide));
    module.querySelectorAll(`svg [data-group="${button.dataset.groupToggle}"]`).forEach(group => group.classList.toggle('hidden-group', hide));
  }));
  });
  document.getElementById('print-report').addEventListener('click', () => window.print());
  document.getElementById('export-summary').addEventListener('click', () => {
    const blob = new Blob([JSON.stringify(result,null,2)],{type:'application/json;charset=utf-8'});
    const url = URL.createObjectURL(blob), link = document.createElement('a');
    link.href = url; link.download = 'Sample01_research_summary_v05.json'; link.click();
    setTimeout(() => URL.revokeObjectURL(url),1000);
    status.textContent = '已导出质控、坐标证据与投影摘要。';
  });
  const progress = document.querySelector('.reading-progress');
  const navLinks = [...document.querySelectorAll('nav a')];
  const navSections = navLinks.map(link=>document.querySelector(link.hash));
  let queued = false;
  function updateProgress() {
    const span = document.documentElement.scrollHeight - innerHeight;
    progress.style.transform = `scaleX(${span > 0 ? Math.max(0,Math.min(1,scrollY/span)) : 0})`;
    let current = navSections[0].id;
    const limit = document.querySelector('.topbar').getBoundingClientRect().height + 120;
    navSections.forEach(section=>{if(section.getBoundingClientRect().top<=limit)current=section.id;});
    navLinks.forEach(link=>{const active=link.hash===`#${current}`;link.classList.toggle('current',active);
      if(active)link.setAttribute('aria-current','location');else link.removeAttribute('aria-current');});
    queued = false;
  }
  addEventListener('scroll', () => {if(!queued){queued=true;requestAnimationFrame(updateProgress);}},{passive:true});
  addEventListener('resize',updateProgress);updateProgress();
  const sections = document.querySelectorAll('section[id]');
  if ('IntersectionObserver' in window) {
    if(!reduce.matches){const observer = new IntersectionObserver(entries=>entries.forEach(entry=>{
      if(entry.isIntersecting){entry.target.classList.add('enter');observer.unobserve(entry.target);}
    }),{threshold:.06});sections.forEach(section=>observer.observe(section));}
  }
})();
