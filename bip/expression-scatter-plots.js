import puppeteer from 'puppeteer'
import { parseArgs } from 'node:util';

const args = process.argv.slice(2);

// const args = ['--gene', 'BRCA1']
const options = {
  gene: {type: 'string'}
};
const { values } = parseArgs({ args, options });
// console.log(values);
// console.log(values.gene);

(async () => {
  const gene = values.gene
  console.log('gene')
  console.log(gene)

  const browser = await puppeteer.launch();
  const page = await browser.newPage();
  await page.setViewport({
    width: 1440,
    height: 1000,
    deviceScaleFactor: 1,
  });
  await page.goto('https://singlecell-staging.broadinstitute.org/single_cell/study/SCP303/human-milk-differential-expression#study-visualize');
  // await page.waitForSelector('svg.gene-load-spinner', {hidden: true})
  console.log('gene again')
  console.log(gene)
  // await page.$eval('.gene-keyword-search input', (el, gene) => el.value = 'BRCA1', gene);
  await page.type('.gene-keyword-search input', gene, {delay: 20})
  await page.keyboard.press('Enter');
  await page.$eval('.gene-keyword-search button', el => el.click());
  await page.waitForSelector('svg.gene-load-spinner', {hidden: true})
  await page.screenshot({path: gene + '.png'});

  await browser.close();
})();
