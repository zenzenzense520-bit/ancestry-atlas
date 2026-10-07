# Design review

The report follows a paper-and-ink editorial direction: restrained cinnabar accents, CJK serif headings, generous column margins, quiet rules, and a specimen index. These choices keep the scientific content legible while giving the report an original visual identity.

| Criterion | Implementation and review |
| --- | --- |
| Typography | Renamed, embedded CJK variable-font subsets; serif headings and tabular numerals; readable mobile body text |
| Whitespace | Consistent section rhythm, wide desktop margins, compact mobile gutters, no horizontal page overflow |
| Hierarchy | Main finding first, separate global/East Asian figures, diagnostics beside each comparison, detailed audits in disclosures |
| Color | Paper, ink and cinnabar; muted reference colors with written group labels and a distinct sample diamond |
| Motion | Short heading entrances, reading progress and state transitions; stationary scientific figures and cards; reduced-motion support |
| Microinteraction | Independent view/group controls, keyboard Space and arrow keys, visible focus, live captions, JSON export feedback and print |
| Responsive layout | Recomputed square scientific figures on narrow screens, fluid desktop figures, one-column cards and contained table scrolling |
| Originality | A sample-specific editorial system built around calculated evidence; open-source projects inform interaction and chart thinking |

Automated checks cover 360, 390, 768 and 1440 px. They verify every computed reference point, unique SVG IDs, independent controls, export, reduced motion, print restoration, absence of background requests and runtime errors. Desktop and mobile hero/figure screenshots are visually reviewed. The report makes no claim to have won an Awwwards, Webby Awards or FWA award.

Reference projects: Observable Plot for layered scientific displays; Lenis for controlled reading rhythm. The implementation uses native JavaScript and Matplotlib SVGs, with no external runtime dependency in the generated HTML.
