"""Offline checks for the archived AlphaFold job builder and ESMFold request."""

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import original_helpers as helpers


class OriginalHelpersTest(unittest.TestCase):
    def test_alphafold_job(self):
        placeholder = "X"
        self.assertEqual(helpers.job("Example-File", placeholder, 1), {
            "name": "example_file",
            "modelSeeds": [],
            "sequences": [{"proteinChain": {"sequence": placeholder, "count": 1}}],
        })

    def test_esmfold_request_and_cached_reply(self):
        reply = "HEADER    OFFLINE TEST\n" + " " * 1100
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mock.pdb"
            with patch.object(helpers.urllib.request, "urlopen",
                              return_value=io.BytesIO(reply.encode())) as mock_open:
                # Close the file explicitly around the archived function's single write.
                with path.open("w") as output, patch.object(helpers, "open", return_value=output, create=True):
                    self.assertTrue(helpers.esmfold("X", str(path)))
                request = mock_open.call_args.args[0]
                self.assertEqual(request.get_method(), "POST")
                self.assertEqual(request.data, b"X")
                self.assertEqual(mock_open.call_args.kwargs["timeout"], 180)
                self.assertEqual(path.read_text(), reply)
                self.assertTrue(helpers.esmfold("X", str(path)))
                mock_open.assert_called_once()


if __name__ == "__main__":
    unittest.main()
