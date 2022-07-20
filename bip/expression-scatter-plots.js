import puppeteer from 'puppeteer'
import { parseArgs } from 'node:util';

const args = process.argv.slice(2);

const options = {
  accession: {type: 'string'}
};
const { values } = parseArgs({ args, options });

(async () => {

  const accession = values.accession
  console.log(`Accession: ${accession}`)

  const browser = await puppeteer.launch();
  const page = await browser.newPage();
  await page.setViewport({
    width: 1680,
    height: 1000,
    deviceScaleFactor: 1,
  });
  const origin = 'https://singlecell-staging.broadinstitute.org'
  const exploreViewUrl = `${origin}/single_cell/study/${accession}#study-visualize`
  const exploreApiUrl =  `${origin}/single_cell/api/v1/studies/${accession}/explore`

  const [exploreApiResponse] = await Promise.all([
    page.waitForResponse(response => response.url() === exploreApiUrl),

    page.goto(exploreViewUrl)
  ])
  // await page.waitForSelector('svg.gene-load-spinner', {hidden: true})

  const exploreJson = await exploreApiResponse.json();

  // All genes in this study
  const uniqueGenes = exploreJson.uniqueGenes

  console.log(`Number of genes: ${uniqueGenes.length}`)

  // Pick a random gene
  const geneIndex = Math.floor(Math.random() * uniqueGenes.length)
  const gene = uniqueGenes[geneIndex]

  // Trigger a gene search
  await page.type('.gene-keyword-search input', gene, {delay: 1})
  await page.keyboard.press('Enter');
  await page.$eval('.gene-keyword-search button', el => el.click());
  console.log(`Waiting for expression plot for gene: ${gene}`)
  const expressionPlotStartTime = Date.now()

  // Wait for UI element signaling that expression plot has finished rendering
  await page.waitForSelector('.cluster-title .badge')
  const expressionPlotPerfTime = Date.now() - expressionPlotStartTime
  console.log(`Expression plot time: ${expressionPlotPerfTime} ms`)

  // Height and width of plot, x- and y-offset from viewport origin
  const clipDimensions = {height: 595, width: 660, x: 5, y: 375}

  // Take a screenshot, save it locally.
  const imagePath = `images/${gene}.webp`
  await page.screenshot({path: imagePath, type: 'webp', clip: clipDimensions});

  console.log(`Wrote ${imagePath}`)

  await browser.close();
})();
