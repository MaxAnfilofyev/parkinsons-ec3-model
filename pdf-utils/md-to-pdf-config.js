
// A marked renderer for mermaid diagrams
const renderer = {
    code(code, infostring) {
        if (infostring === 'mermaid'){
            return `<pre class="mermaid">${code}</pre>`
        }
        return false
    },
};

module.exports = {
    marked_extensions: [{ renderer }],
    script: [
		{ url: 'https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js' },
		{ content: 'mermaid.initialize({ startOnLoad: false}); (async () => { await mermaid.run(); })();' },
		{url: 'http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'},
		{content: 'MathJax.Hub.Config({ extensions: ["tex2jax.js"], jax: ["input/TeX", "output/HTML-CSS"], tex2jax: { inlineMath: [ ["$","$"], ["\\\\(","\\\\)"] ], displayMath: [ ["$$","$$"], ["\\\\[","\\\\]"] ], }, "HTML-CSS": { availableFonts: ["TeX"] } });'}

	],
	pdf_options: {
		format: 'Letter',
		margin: '20mm',
		printBackground: true,
		headerTemplate: 
			'<div style="font-size:9px;width:100%;text-align:center;margin-bottom:5px;"><span class="date"></span></div>',
		footerTemplate: 
			'<div style="font-size:9px;width:100%;text-align:left;margin-bottom:5px;margin-left:20mm"><span class="pageNumber"></span> of <span class="totalPages"></span>   -   bioRχiv</div><div style="font-size:9px;width:100%;text-align:right;margin-bottom:5px;margin-right:20mm">Anfilofyev   -   Energetic Tipping Point in SNc Neurons</div>',
		
	}
};
