const fs = require('fs');
const spawn = require('child_process').spawnSync;

//SAINT file to process
const matrix = process.argv[2];
const rankStart = process.argv[3];
const rankTotal = process.argv[4];

// parameters
const resultsFolder = 'results';

// create Results folder if it doesn't exists
if (!fs.existsSync(resultsFolder)) {
	fs.mkdirSync(resultsFolder);
}

// ranks to test
const ranks = [...Array(rankTotal).keys()].map(x => x + rankStart);

// NMF
const nmf = function() {
	ranks.forEach((rank) => {
		const outFolder = `Results/${rank}rank`;
		fs.mkdirSync(outFolder);
		const nmfProcess = spawn('python3', [
			`${__dirname}/nmf/nmf_bcv/nmf.py`,
			'--input', matrix,
			'--k', rank,
			'--model-output', `${outFolder}/model.pkl`,
			'--basis-output', `${outFolder}/basis.csv`,
			'--score-output', `${outFolder}/scores.csv`,
			'--l1-ratio', 1
		]);
		console.log(`completed NMF for rank ${rank}`);
		nmfSummary(outFolder);
	});
}

// NMF sumarization
const nmfSummary = function(outFolder) {
	const rProcess = spawn(`${__dirname}/nmf/rankProfile.R`, [
		`${outFolder}/basis.csv`,
		`${outFolder}/scores.csv`,
		'GO:CC',
		outFolder
	]);
	console.log(`completed NMF summarization`);
}

nmf();
