module.exports = {
    "env": {
        "es6": true,
        "jest": true,
        "node": true
    },
    "extends": "eslint:recommended",
    "parserOptions": {
        "ecmaVersion": 2016,
        "ecmaFeatures": {
          "experimentalObjectRestSpread": true
        }
    },
    "rules": {
        "comma-dangle": ["error", {
          "arrays": "always-multiline",
          "objects": "always-multiline",
          "imports": "always-multiline",
          "exports": "always-multiline",
          "functions": "ignore"
        }],
        "eol-last": [
          "error",
          "always"
        ],
        "indent": [
            "error",
            2
        ],
        "linebreak-style": [
            "error",
            "unix"
        ],
        "max-len": [
          "error",
          { "code": 100 },
        ],
        "no-trailing-spaces": ["error"],
        "quotes": [
            "error",
            "single"
        ],
        "semi": [
            "error",
            "always"
        ]
    }
};